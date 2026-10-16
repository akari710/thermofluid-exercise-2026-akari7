# 実験・診断・保存の配布実装。includeだけでは実行しない。
using Plots

function integrate(f,y; method="forward_euler",N=8,T=1.0,A=nothing,
                   factorizer=lu,linear_solver=(F,b)->F\b)
    method in method_ids(true) || throw(ArgumentError("未知のmethod ID: $method"))
    N isa Integer && !(N isa Bool) && N>0 || throw(ArgumentError("Nは正整数が必要です"))
    finite_real(T) && T>0 || throw(ArgumentError("終端Tは有限正が必要です"))
    dt = T/N
    validate_step(y,zero(T),dt)
    implicit = method in ("backward_euler","crank_nicolson")
    implicit && validate_linear(A,y)
    counts = Dict("rhs_evaluations"=>0,"factorizations"=>0,"linear_solves"=>0)
    counted_f(u,t) = (counts["rhs_evaluations"]+=1; f(u,t))
    # 失敗の文脈を保持し、成功結果へ変換しない。
    context(step,t,operation) = try
        operation()
    catch exception
        error("method=$method step=$step time=$t dt=$dt: $(sprint(showerror,exception))")
    end
    stepper = if implicit
        context(0,0.0,()->prepare_linear_step(A,y,0.0,dt;
            theta=method=="backward_euler" ? 1.0 : 0.5,
            factorizer=M->(counts["factorizations"]+=1;factorizer(M)),
            linear_solver=(F,b)->(counts["linear_solves"]+=1;linear_solver(F,b))))
    else
        step = method=="forward_euler" ? forward_euler_step : method=="midpoint_rk2" ? rk2_step : classical_rk4_step
        (u,t)->step(counted_f,u,t,dt)
    end
    states = y isa Real ? Real[state_copy(y)] : [state_copy(y)]
    times = collect(0:N).*dt
    for n in 1:N
        next = context(n,times[n],()->begin
            value = stepper(states[end],times[n])
            valid_state(value) || error("更新結果が非有限値です")
            value
        end)
        push!(states,state_copy(next))
    end
    (;states,times,dt,N,final_time=times[end],cost=counts)
end

function experiments(;include_cn=false)
    include_cn isa Bool || throw(ArgumentError("include_cnはBoolが必要です"))
    convergence = Dict{String,Any}[]
    for method in method_ids(include_cn)
        runs = [integrate((y,t)->-y,1.0;method,N,T=1.0,A=-1.0) for N in (8,16,32,64)]
        errors = [abs(r.states[end]-exp(-1)) for r in runs]
        all(e->isfinite(e)&&e>0,errors) || error("誤差は有限正が必要です: $method")
        push!(convergence,Dict("method"=>method,"N"=>[r.N for r in runs],
            "dt"=>[r.dt for r in runs],"final_time"=>[r.final_time for r in runs],
            "final_state"=>[r.states[end] for r in runs],"errors"=>errors,
            "orders"=>log2.(errors[1:end-1]./errors[2:end]),"cost"=>[r.cost for r in runs]))
    end
    grid = diffusion_problem()
    diffusion = Dict{String,Any}[]
    for N in (128,4), method in method_ids(include_cn)
        r = integrate((y,t)->grid.A*y,grid.initial;method,N,T=0.1,A=grid.A)
        reference = exp(grid.lambda1*r.final_time)*grid.phi1 + 0.01exp(grid.lambda39*r.final_time)*grid.phi39
        high = [dot(u,grid.phi39)/dot(grid.phi39,grid.phi39) for u in r.states]
        push!(diffusion,Dict("method"=>method,"N"=>N,"dt"=>r.dt,"final_time"=>r.final_time,
            "Fo"=>grid.D*r.dt/grid.dx^2,"times"=>r.times,"high_mode_coefficients"=>high,
            "initial_state"=>vcat(0.0,grid.initial,0.0),"final_state"=>vcat(0.0,r.states[end],0.0),
            "reference_state"=>vcat(0.0,reference,0.0),
            "max_abs_error"=>maximum(abs.(r.states[end]-reference)),"cost"=>r.cost))
    end
    Dict("schema_version"=>1,"example_id"=>"time-integration",
        "time_problem"=>Dict("rhs"=>"-y","initial_state"=>1.0,"T"=>1.0,"reference"=>exp(-1)),
        "time_convergence"=>convergence,
        "diffusion_problem"=>Dict("domain"=>[0.0,2.0],"D"=>grid.D,"dx"=>grid.dx,"x"=>grid.x,
            "points"=>41,"interior_unknowns"=>39,"boundary"=>"Dirichlet 0 at both ends",
            "space_scheme"=>"second-order centered difference","T"=>0.1,
            "initial_modes"=>[1,39],"initial_coefficients"=>[1.0,0.01],
            "eigenvalues"=>[grid.lambda1,grid.lambda39],"reference"=>"semidiscrete exact solution"),
        "diffusion"=>diffusion)
end

const OUTPUT_FILES = ("time-convergence.png","diffusion-comparison.png","summary.toml")
default_output_dir() = normpath(joinpath(@__DIR__,"..","..","scratch","time-integration"))

function validate_output_dir(output_dir)
    output_dir isa AbstractString && !isempty(strip(output_dir)) ||
        throw(ArgumentError("出力先は空でないパスが必要です"))
    directory = abspath(output_dir)
    (ispath(directory) || islink(directory)) && (!isdir(directory) || islink(directory)) &&
        throw(ArgumentError("出力先は通常のディレクトリが必要です"))
    # 書込み前に既存ファイルの形を検査する。無関係なディレクトリを置換しない。
    parent = dirname(directory)
    while !ispath(parent) && !islink(parent)
        parent = dirname(parent)
    end
    isdir(parent) || throw(ArgumentError("出力先の親がディレクトリではありません"))
    for name in OUTPUT_FILES
        path = joinpath(directory,name)
        (ispath(path) || islink(path)) && (!isfile(path) || islink(path)) &&
            throw(ArgumentError("出力ファイルは通常ファイルが必要です: $path"))
    end
    directory
end

function render_plots(data,directory)
    ids = method_ids(true)
    colors = (:dodgerblue3,:darkorange2,:forestgreen,:purple3,:brown3)
    styles = (:solid,:dash,:dot,:dashdot,:dashdotdot)
    color(id) = colors[findfirst(==(id),ids)]
    style(id) = styles[findfirst(==(id),ids)]
    p = plot(;xscale=:log10,yscale=:log10,xlabel="dt",ylabel="Absolute error at T = 1",
        xticks=([1/64,1/32,1/16,1/8],["1/64","1/32","1/16","1/8"]),
        title="Time convergence: y' = -y",legend=:topleft,size=(900,600))
    for row in data["time_convergence"]
        id = row["method"]
        plot!(p,row["dt"],row["errors"];marker=:circle,color=color(id),linestyle=style(id),label=id)
    end
    savefig(p,joinpath(directory,OUTPUT_FILES[1]))
    x = data["diffusion_problem"]["x"]
    lambda39 = data["diffusion_problem"]["eigenvalues"][2]
    panels = []
    for N in (128,4)
        rows = filter(r->r["N"]==N,data["diffusion"])
        fo = first(rows)["Fo"]
        final = plot(x,first(rows)["reference_state"];color=:black,linestyle=:dash,linewidth=2,
            label="semidiscrete exact",xlabel="x",ylabel="u(x, 0.1)",title="Final: N=$N, Fo=$fo",legend=:topright)
        mode = plot(;xlabel="t",ylabel="Mode 39 coefficient",title="All methods: N=$N",
            yscale=N==128 ? :log10 : :identity,legend=:topleft)
        times = first(rows)["times"]
        implicit = plot(times,0.01exp.(lambda39.*times);color=:black,linestyle=:dash,
            label="semidiscrete exact",xlabel="t",ylabel="Mode 39 coefficient",
            title="Implicit methods: N=$N",yscale=N==128 ? :log10 : :identity,legend=:topright)
        for row in rows
            id = row["method"]
            plot!(final,x,row["final_state"];label=id,color=color(id),linestyle=style(id))
            plot!(mode,row["times"],row["high_mode_coefficients"];label=id,color=color(id),
                linestyle=style(id),marker=N==4 ? :circle : :none)
            if id in ("backward_euler","crank_nicolson")
                plot!(implicit,row["times"],row["high_mode_coefficients"];label=id,color=color(id),
                    linestyle=style(id),marker=N==4 ? :circle : :none)
            end
        end
        # 陰解法だけの表示で、小さいCNの符号交代を陽解法の増幅と同時に読める。
        push!(panels,final,mode,implicit)
    end
    p = plot(panels...;layout=(2,3),size=(1800,1000),
        bottom_margin=6Plots.mm,left_margin=12Plots.mm,legendfontsize=7)
    savefig(p,joinpath(directory,OUTPUT_FILES[2]))
end

write_summary(data,path) = open(io->TOML.print(io,data;sorted=true),path,"w")

# 3出力を一時領域で完成させる。反映失敗時は旧ファイルを復元する。
# hooksは教員側の描画・保存・反映失敗検査にも使う。
function save_outputs(data,output_dir;plot_writer=render_plots,summary_writer=write_summary,
                      install_file=(source,target)->mv(source,target;force=true))
    directory = validate_output_dir(output_dir)
    mktempdir() do stage
        plot_writer(data,stage)
        summary_writer(data,joinpath(stage,OUTPUT_FILES[3]))
        all(name->isfile(joinpath(stage,name)) && filesize(joinpath(stage,name))>0,OUTPUT_FILES) ||
            error("3出力が完成していません")
        backup = joinpath(stage,"previous"); mkpath(backup)
        existed = Dict(name=>isfile(joinpath(directory,name)) for name in OUTPUT_FILES)
        for name in OUTPUT_FILES
            existed[name] && cp(joinpath(directory,name),joinpath(backup,name))
        end
        mkpath(directory)
        try
            for name in OUTPUT_FILES
                install_file(joinpath(stage,name),joinpath(directory,name))
            end
        catch
            for name in OUTPUT_FILES
                target = joinpath(directory,name)
                if existed[name]
                    cp(joinpath(backup,name),target;force=true)
                elseif isfile(target)
                    rm(target)
                end
            end
            rethrow()
        end
    end
    directory
end

function main(;output_dir=default_output_dir(),include_cn=false,experiment_runner=experiments,kwargs...)
    directory = validate_output_dir(output_dir)
    include_cn isa Bool || throw(ArgumentError("include_cnはBoolが必要です"))
    data = experiment_runner(;include_cn)
    save_outputs(data,directory;kwargs...)
    data
end
