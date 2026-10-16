# 実験・診断・描画・保存の配布実装。学生の必須編集箇所ではありません。
# run.jlから読み込み、TimeIntegrationExampleモジュールの一部になります。
# 読む順: 共通の時間ループ → 実験A/B → 出力先の検証 → 描画 → 保存 → main。
# includeだけでは実験・出力を行いません。
using Plots

"""
    integrate(f, y; method="forward_euler", N=8, T=1.0, A=nothing,
              factorizer=lu, linear_solver=...)

共通の時間ループで、初期状態 `y` を時刻0から `T` まで `N` ステップ進める。

# 引数

- `f`: 陽解法で使う右辺 `f(y, t)`。
- `y`: 有限実数スカラー、または非空で全要素が有限な `Vector{Float64}`。
- `method`: `method_ids` にある固定ID。
- `N`, `T`: 正整数のステップ数と、有限な正の終端時刻。
- `A`: 陰解法で使う線形問題の係数。陽解法では使わない。
- `factorizer`, `linear_solver`: 分解・求解の関数。通常は既定値を使う。

入力は計算前に検証し、時間刻みは `dt = T/N` に固定する。
陽解法の安定性のためにdtを自動縮小しない。
計算失敗時は、手法・ステップ・時刻・dtを示して停止する。

# 返り値

全ステップの `states` と `times`、`dt`、`N`、実終端時刻 `final_time`、
積分中の実際のRHS評価・分解・求解回数 `cost` を持つ名前付きタプル。
入力状態は変更せず、各ステップの状態を独立して保存する。
"""
function integrate(f,y; method="forward_euler",N=8,T=1.0,A=nothing,
                   factorizer=lu,linear_solver=(F,b)->F\b)
    # 時間ループへ入る前に、手法・ステップ数・終端と初期状態を検証する。
    method in method_ids(true) || throw(ArgumentError("未知のmethod ID: $method"))
    N isa Integer && !(N isa Bool) && N>0 || throw(ArgumentError("Nは正整数が必要です"))
    finite_real(T) && T>0 || throw(ArgumentError("終端Tは有限正が必要です"))

    # 全手法の終端を揃える。陰解法では状態と係数Aの対応も調べる。
    dt = T/N
    validate_step(y,zero(T),dt)
    implicit = method in ("backward_euler","crank_nicolson")
    implicit && validate_linear(A,y)

    # 呼出し時に回数を増やし、段数から推測した値で代用しない。
    counts = Dict("rhs_evaluations"=>0,"factorizations"=>0,"linear_solves"=>0)
    counted_f(u,t) = (counts["rhs_evaluations"]+=1; f(u,t))

    # 失敗の文脈を保持し、成功結果へ変換しない。
    context(step,t,operation) = try
        operation()
    catch exception
        error("method=$method step=$step time=$t dt=$dt: $(sprint(showerror,exception))")
    end

    # 陰解法は時間ループの外で1回だけ分解し、返された関数を使い回す。
    stepper = if implicit
        context(0,0.0,()->prepare_linear_step(A,y,0.0,dt;
            theta=method=="backward_euler" ? 1.0 : 0.5,
            factorizer=M->(counts["factorizations"]+=1;factorizer(M)),
            linear_solver=(F,b)->(counts["linear_solves"]+=1;linear_solver(F,b))))
    else
        # 陽解法は同じRHSとdtを使い、選んだ手法の1ステップ関数を呼ぶ。
        step = method=="forward_euler" ? forward_euler_step : method=="midpoint_rk2" ? rk2_step : classical_rk4_step
        (u,t)->step(counted_f,u,t,dt)
    end

    # 初期状態を含むN+1時刻を保存する。times[n]はn回目の更新前の時刻。
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

"""
    experiments(; include_cn=false)

実験Aの時間収束と、実験Bの拡散安定性を計算し、結果の辞書を返す。
`include_cn=true` なら任意のCrank–Nicolsonも比較する。

実験Aは常微分方程式 `y′ = -y`、実験Bは同じ空間離散化の拡散方程式を使う。
誤差・次数・モード係数の診断は、積分中の計算回数に含めない。
この関数はファイルを保存しない。
返り値の構造は `summary.toml` に対応し、格子・初期条件・実際のdtも記録する。
RK2が未実装なら、そのエラーで停止する。
"""
function experiments(;include_cn=false)
    include_cn isa Bool || throw(ArgumentError("include_cnはBoolが必要です"))

    # 実験A: 空間誤差を含まない問題で、時間刻みを半分ずつにする。
    convergence = Dict{String,Any}[]
    for method in method_ids(include_cn)
        runs = [integrate((y,t)->-y,1.0;method,N,T=1.0,A=-1.0) for N in (8,16,32,64)]

        # 同じT=1での絶対誤差から、隣り合う2つの刻みの観測次数を求める。
        errors = [abs(r.states[end]-exp(-1)) for r in runs]
        all(e->isfinite(e)&&e>0,errors) || error("誤差は有限正が必要です: $method")
        push!(convergence,Dict("method"=>method,"N"=>[r.N for r in runs],
            "dt"=>[r.dt for r in runs],"final_time"=>[r.final_time for r in runs],
            "final_state"=>[r.states[end] for r in runs],"errors"=>errors,
            "orders"=>log2.(errors[1:end-1]./errors[2:end]),"cost"=>[r.cost for r in runs]))
    end

    # 実験B: Nごとに格子・初期値・dtを揃え、手法だけを変える。
    grid = diffusion_problem()
    diffusion = Dict{String,Any}[]
    for N in (128,4), method in method_ids(include_cn)
        r = integrate((y,t)->grid.A*y,grid.initial;method,N,T=0.1,A=grid.A)

        # 参照も同じ差分行列の固有値を使い、時間誤差だけを比較する。
        reference = exp(grid.lambda1*r.final_time)*grid.phi1 + 0.01exp(grid.lambda39*r.final_time)*grid.phi39

        # dot(u,phi)/dot(phi,phi)でモード39の係数を取り出す。
        # 各時刻の係数を並べ、高周波の増幅・符号交代を追う。
        high = [dot(u,grid.phi39)/dot(grid.phi39,grid.phi39) for u in r.states]

        # 内部39点の計算結果に、既知の両端0を足して41点の分布に戻す。
        push!(diffusion,Dict("method"=>method,"N"=>N,"dt"=>r.dt,"final_time"=>r.final_time,
            "Fo"=>grid.D*r.dt/grid.dx^2,"times"=>r.times,"high_mode_coefficients"=>high,
            "initial_state"=>vcat(0.0,grid.initial,0.0),"final_state"=>vcat(0.0,r.states[end],0.0),
            "reference_state"=>vcat(0.0,reference,0.0),
            "max_abs_error"=>maximum(abs.(r.states[end]-reference)),"cost"=>r.cost))
    end

    # 数値結果と再現に必要な条件を、描画・TOML保存へ渡す。
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

# この3ファイルを1組として保存する。無関係な出力は扱わない。
const OUTPUT_FILES = ("time-convergence.png","diffusion-comparison.png","summary.toml")

"""
    default_output_dir()

学生リポジトリ内の `scratch/time-integration/` の絶対パスを返す。
現在の作業ディレクトリではなく、このファイルの場所を基準に決める。
"""
default_output_dir() = normpath(joinpath(@__DIR__,"..","..","scratch","time-integration"))

"""
    validate_output_dir(output_dir)

出力先を検証し、絶対パスを返す。ファイル・ディレクトリは作成しない。
空のパス、ディレクトリでない出力先、3出力の位置にあるディレクトリや
シンボリックリンクなどは `ArgumentError` で拒否する。
既存の成果物や無関係なファイルを、入力検証で変更しない。
"""
function validate_output_dir(output_dir)
    output_dir isa AbstractString && !isempty(strip(output_dir)) ||
        throw(ArgumentError("出力先は空でないパスが必要です"))

    # 出力先そのものが、通常ファイルやリンクでないことを確認する。
    directory = abspath(output_dir)
    (ispath(directory) || islink(directory)) && (!isdir(directory) || islink(directory)) &&
        throw(ArgumentError("出力先は通常のディレクトリが必要です"))

    # 書込み前に既存ファイルの形を検査する。無関係なディレクトリを置換しない。
    parent = dirname(directory)
    while !ispath(parent) && !islink(parent)
        parent = dirname(parent)
    end
    isdir(parent) || throw(ArgumentError("出力先の親がディレクトリではありません"))

    # 置換対象3点だけを検査する。他の名前のファイルは保持する。
    for name in OUTPUT_FILES
        path = joinpath(directory,name)
        (ispath(path) || islink(path)) && (!isfile(path) || islink(path)) &&
            throw(ArgumentError("出力ファイルは通常ファイルが必要です: $path"))
    end
    directory
end

"""
    render_plots(data, directory)

`experiments` の計算済み結果からPNG2点を `directory` に保存する。
時間収束は両対数軸、拡散は最終分布・全手法のモード39・陰解法のモード39を描く。
粗いdtのモード係数は線形軸で描き、CNの負値も見えるようにする。
ここでは再積分せず、渡されたデータから図を作る。
"""
function render_plots(data,directory)
    # どのパネルでも、同じ手法には同じ色と線種を使う。
    ids = method_ids(true)
    colors = (:dodgerblue3,:darkorange2,:forestgreen,:purple3,:brown3)
    styles = (:solid,:dash,:dot,:dashdot,:dashdotdot)
    color(id) = colors[findfirst(==(id),ids)]
    style(id) = styles[findfirst(==(id),ids)]
    names = ("Forward Euler","Midpoint RK2","Classical RK4","Backward Euler","Crank-Nicolson")
    name(id) = names[findfirst(==(id),ids)]

    # 実験A: 横軸dt、縦軸誤差。直線の傾きが時間積分の次数に対応する。
    p = plot(;xscale=:log10,yscale=:log10,xlabel="dt",ylabel="Absolute error at T = 1",
        xticks=([1/64,1/32,1/16,1/8],["1/64","1/32","1/16","1/8"]),
        title="Time convergence: y' = -y",legend=:bottomright,size=(900,600),
        left_margin=8Plots.mm,bottom_margin=6Plots.mm,legendfontsize=11)
    for row in data["time_convergence"]
        id = row["method"]
        plot!(p,row["dt"],row["errors"];marker=:circle,color=color(id),linestyle=style(id),label=name(id))
    end
    savefig(p,joinpath(directory,OUTPUT_FILES[1]))

    # 実験B: 左列を細かいdt、右列を粗いdtとし、比較する量を上下に揃える。
    x = data["diffusion_problem"]["x"]
    lambda39 = data["diffusion_problem"]["eigenvalues"][2]
    markers = (:circle,:rect,:diamond,:utriangle,:dtriangle)
    marker(id) = markers[findfirst(==(id),ids)]
    columns = []
    for (column,N) in enumerate((128,4))
        rows = filter(r->r["N"]==N,data["diffusion"])
        fo = first(rows)["Fo"]
        T = first(rows)["final_time"]
        time_limits = (-0.03*T,1.03*T)
        reference_times = range(0.0,T;length=401)
        reference_mode = 0.01exp.(lambda39.*reference_times)
        letters = column==1 ? ("a","c","e") : ("b","d","f")

        # 上から最終分布、全手法の高周波係数、陰解法の高周波係数。
        final = plot(;xlabel="x",ylabel="u(x, 0.1)",
            title="($(letters[1])) Final profile\nN = $N, Fo = $fo",legend=false)
        mode = plot(;xlabel="t",ylabel="Mode 39 coefficient",
            title="($(letters[2])) All methods: N = $N",
            yscale=N==128 ? :log10 : :identity,legend=false)
        implicit = plot(;xlabel="t",ylabel="Mode 39 coefficient",
            title="($(letters[3])) Implicit methods: N = $N",
            yscale=N==128 ? :log10 : :identity,legend=false)
        for panel in (mode,implicit)
            plot!(panel;xlims=time_limits,
                xticks=([0.0,0.025,0.05,0.075,0.1],["0","0.025","0.05","0.075","0.1"]))
            N==4 && hline!(panel,[0.0];color=:gray65,linewidth=1,label="")
        end

        for row in rows
            id = row["method"]
            plot!(final,x,row["final_state"];label="",color=color(id),linestyle=style(id),linewidth=2)
            plot!(mode,row["times"],row["high_mode_coefficients"];label="",color=color(id),
                linestyle=style(id),linewidth=2.2,marker=N==4 ? marker(id) : :none,
                markersize=5,markercolor=:white,markerstrokecolor=color(id),markerstrokewidth=1.5)
            if id in ("backward_euler","crank_nicolson")
                plot!(implicit,row["times"],row["high_mode_coefficients"];label="",color=color(id),
                    linestyle=style(id),linewidth=2.2,marker=N==4 ? marker(id) : :none,
                    markersize=5,markercolor=:white,markerstrokecolor=color(id),markerstrokewidth=1.5)
            end
        end

        plot!(final,x,first(rows)["reference_state"];color=:black,linestyle=:dash,linewidth=2,label="")
        for panel in (mode,implicit)
            plot!(panel,reference_times,reference_mode;color=:black,linestyle=:dash,linewidth=2,label="")
        end
        push!(columns,(final,mode,implicit))
    end

    # 共通凡例を上にまとめ、長い線見本で線種を判別できるようにする。
    key = plot(;xlims=(0,3),ylims=(0,2),showaxis=false,grid=false,legend=false,
        ticks=false,top_margin=0Plots.mm,bottom_margin=0Plots.mm)
    active_ids = [row["method"] for row in data["time_convergence"]]
    for (i,id) in enumerate(vcat(["reference"],active_ids))
        col = (i-1)%3
        y = 1.5-div(i-1,3)
        exact = id=="reference"
        plot!(key,[col+0.02,col+0.24],[y,y];color=exact ? :black : color(id),
            linestyle=exact ? :dash : style(id),linewidth=2.2,label="")
        annotate!(key,col+0.29,y,text(exact ? "Semidiscrete exact" : name(id),11,:left))
    end
    panels = [columns[col][row] for row in 1:3 for col in 1:2]
    arrangement = @layout [key{0.08h}; grid(3,2)]
    p = plot(key,panels...;layout=arrangement,size=(1200,1500),
        bottom_margin=6Plots.mm,left_margin=8Plots.mm,right_margin=4Plots.mm,
        titlefontsize=13,guidefontsize=12,tickfontsize=10,gridalpha=0.15)
    savefig(p,joinpath(directory,OUTPUT_FILES[2]))
end

"""
    write_summary(data, path)

計算結果の辞書をTOML形式で `path` に書き込む。
通常は `save_outputs` が用意する一時ファイルに対して呼び出す。
"""
write_summary(data,path) = open(io->TOML.print(io,data;sorted=true),path,"w")

"""
    save_outputs(data, output_dir; plot_writer=render_plots,
                 summary_writer=write_summary, install_file=...)

PNG2点とTOML1点を一時領域で完成させ、1組の成果物として出力先へ反映する。

# 引数

- `data`: `experiments` の計算結果。
- `output_dir`: 出力先。省略はできず、保存前に検証する。
- `plot_writer`, `summary_writer`, `install_file`: 作図・TOML保存・反映の関数。
  通常は既定値を使い、教員の検査では失敗を起こす関数へ差し替える。

作図・TOML保存に失敗した場合は既存出力を変更しない。
反映途中で失敗した場合は、退避した旧出力を復元して同じ失敗を再通知する。
3点以外のファイルは変更しない。

# 返り値

反映に成功した出力先の絶対パス。
"""
function save_outputs(data,output_dir;plot_writer=render_plots,summary_writer=write_summary,
                      install_file=(source,target)->mv(source,target;force=true))
    directory = validate_output_dir(output_dir)

    # 作図とTOML保存が両方成功するまで、出力先へ書き込まない。
    mktempdir() do stage
        plot_writer(data,stage)
        summary_writer(data,joinpath(stage,OUTPUT_FILES[3]))
        all(name->isfile(joinpath(stage,name)) && filesize(joinpath(stage,name))>0,OUTPUT_FILES) ||
            error("3出力が完成していません")

        # 反映直前に旧3点を退避し、もともとなかった出力も区別する。
        backup = joinpath(stage,"previous"); mkpath(backup)
        existed = Dict(name=>isfile(joinpath(directory,name)) for name in OUTPUT_FILES)
        for name in OUTPUT_FILES
            existed[name] && cp(joinpath(directory,name),joinpath(backup,name))
        end

        # 完成した3点を反映する。途中失敗したら、下のcatchで元へ戻す。
        mkpath(directory)
        try
            for name in OUTPUT_FILES
                install_file(joinpath(stage,name),joinpath(directory,name))
            end
        catch
            # 旧ファイルは復元し、新しくできたファイルだけなら取り除く。
            for name in OUTPUT_FILES
                target = joinpath(directory,name)
                if existed[name]
                    cp(joinpath(backup,name),target;force=true)
                elseif isfile(target)
                    rm(target)
                end
            end

            # 復元できても処理の成功とは扱わず、元の失敗を呼出し元へ伝える。
            rethrow()
        end
    end
    directory
end

"""
    main(; output_dir=default_output_dir(), include_cn=false,
         experiment_runner=experiments, kwargs...)

出力先を検証し、実験A/Bの計算と3点の保存を順に行う。
`include_cn=true` なら任意CNを追加し、`output_dir` で保存先を指定できる。
`experiment_runner` と `kwargs...` は教員の失敗時検査にも使う。
通常の実行では既定値のままとする。

計算が失敗した場合は保存へ進まず、既存出力を保持する。
返り値は計算結果の辞書。course_progress.tomlは更新しない。
"""
function main(;output_dir=default_output_dir(),include_cn=false,experiment_runner=experiments,kwargs...)
    directory = validate_output_dir(output_dir)
    include_cn isa Bool || throw(ArgumentError("include_cnはBoolが必要です"))

    # 計算をすべて終えてから保存へ進む。未実装RK2もここで停止する。
    data = experiment_runner(;include_cn)
    save_outputs(data,directory;kwargs...)

    data
end
