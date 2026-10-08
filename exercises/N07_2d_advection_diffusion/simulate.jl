# N07: セル中心格子 → 保存区間の更新 → 温度とBurgersの同一実行でのHDF5保存。

module N07Simulation
using ThermofluidExercise
include("provided_support.jl")
const N = ThermofluidExercise.N07

"""
    simulate_case(
        id,
        nx,
        ny;
        times = SAVE_TIMES,
        thermal_step = N.thermal_step!,
        burgers_step = N.burgers_step!,
        burgers_timestep = N.burgers_stable_timestep,
    )

各保存区間で温度または二成分速度を更新し、場と熱輸送を記録する。

# 引数

- `id`: 格子点数と対応するケース名。
- `nx`: x方向の格子点数。
- `ny`: y方向の格子点数。
- `times`: 先頭0、末尾t_finalの有限な狭義増加の保存時刻列。
- `thermal_step`: 新温度と辺別熱輸送レートを返す更新関数。
- `burgers_step`: 旧二成分速度を保持して新二成分を更新する関数。
- `burgers_timestep`: 現在の旧二成分場から刻み上限を返す関数。

# 返り値

条件・座標・保存時刻・全step時刻/刻み・保存step対応と場の辞書。温度は辺別区間積分も含む。時間が進まない刻みや非有限値は条件と時刻を示して停止する。
"""
function simulate_case(
    id,
    nx,
    ny;
    times = SAVE_TIMES,
    thermal_step = N.thermal_step!,
    burgers_step = N.burgers_step!,
    burgers_timestep = N.burgers_stable_timestep,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        nx isa Integer &&
            ny isa Integer &&
            !(nx isa Bool || ny isa Bool) &&
            nx >= 3 &&
            ny >= 3,
        "各軸3以上の整数が必要です",
    )
    time = schedule(times)
    c = case_metadata(id, nx, ny)
    dx, dy = c["dx"], c["dy"]

    # 座標と格子幅を用意し、配列の添字との対応をそろえる。
    x = cell_centers(nx, 2.)
    y = cell_centers(ny, 1.)
    isburgers = id == "periodic_burgers"
    a, v = isburgers ? exact_burgers(x, y, 0.) : (exact_temperature(x, y, 0., id), nothing)
    b = similar(a)
    w = isburgers ? similar(v) : nothing
    nt = length(time)

    # 初期値を保持し、更新に使う作業用配列を用意する。
    U = Array{Float64}(undef, nx, ny, nt)
    U[:, :, 1] = a
    V = isburgers ? similar(U) : nothing
    isburgers && (V[:, :, 1] = v)
    adv = zeros(4, nt - 1)
    diff = zeros(4, nt - 1)
    step_time = Float64[]
    step_dt = Float64[]
    indices = [0]
    t = 0.
    config = isburgers ? nothing : thermal_config(id)

    # 安定上限と到達する時刻に合わせて時間刻みを決める。
    dtmax =
        isburgers ? NaN :
        N.thermal_stable_timestep(config.cx, config.cy, config.kappa, dx, dy, config.bc)
    if startswith(id, "comparison_")
        ref = thermal_config("periodic_combined")
        dtmax = N.thermal_stable_timestep(ref.cx, ref.cy, ref.kappa, dx, dy, ref.bc)
    end

    # 旧場から更新し、診断・境界処理を終えてから次の反復へ進む。
    for k in 2:nt
        target = time[k]
        while t < target
            dt = min(isburgers ? burgers_timestep(a, v, 0.05, dx, dy) : dtmax, target - t)
            require(
                dt isa Real && isfinite(dt) && dt > 0 && t + dt > t,
                "$id t=$t: 時間が進まないdtです",
            )
            try
                if isburgers
                    burgers_step(b, w, a, v, dt, dx, dy; nu = 0.05)
                    require(all(isfinite, b) && all(isfinite, w), "非有限速度です")
                    v, w = w, v
                else
                    r = thermal_step(
                        b,
                        a,
                        dt,
                        dx,
                        dy;
                        cx = config.cx,
                        cy = config.cy,
                        kappa = config.kappa,
                        bc = config.bc,
                    )
                    require(all(isfinite, b), "非有限温度です")
                    for (j, side) in enumerate(SIDES)
                        adv[j, k - 1]+=dt *
                                       getproperty(r.boundary_rates.advective, Symbol(side))
                        diff[j, k - 1]+=dt * getproperty(
                            r.boundary_rates.diffusive,
                            Symbol(side),
                        )
                    end
                end
            catch e
                error("$id $(nx)x$(ny) t=$t dt=$dt: $(sprint(showerror,e))")
            end
            push!(step_time, t)
            push!(step_dt, dt)
            t = min(target, t + dt)
            a, b = b, a
        end
        U[:, :, k] = a
        isburgers && (V[:, :, k] = v)
        push!(indices, length(step_dt))
    end
    merge!(
        c,
        Dict(
            "x" => x,
            "y" => y,
            "time" => time,
            "step_time" => step_time,
            "step_dt" => step_dt,
            "save_step_index" => indices,
        ),
    )
    if isburgers
        c["u"] = U
        c["v"] = V
    else
        c["temperature"] = U
        c["advective_heat_integral"] = adv
        c["diffusive_heat_integral"] = diff
    end
    c
end

"""
    main(;
        output_dir = DEFAULT_OUTPUT_DIR,
        thermal_step = N.thermal_step!,
        burgers_step = N.burgers_step!,
        burgers_timestep = N.burgers_stable_timestep,
        publish_options...,
    )

公式ケースの温度とBurgersを同一run_idで保存する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `thermal_step`: 新温度と辺別熱輸送レートを返す更新関数。
- `burgers_step`: 旧二成分速度を保持して新二成分を更新する関数。
- `burgers_timestep`: 現在の旧二成分場から刻み上限を返す関数。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

nothing。temperature.h5とburgers.h5を一時生成・相互検証してから復元付きで同時反映する。
"""
function main(;
    output_dir = DEFAULT_OUTPUT_DIR,
    thermal_step = N.thermal_step!,
    burgers_step = N.burgers_step!,
    burgers_timestep = N.burgers_stable_timestep,
    publish_options...,
)
    run_id = bytes2hex(rand(UInt8, 16))

    # 一時領域で処理を完了してから、公式出力を反映する。
    staged(output_dir, ("temperature.h5", "burgers.h5"); publish_options...) do stage
        for family in ("temperature", "burgers")
            cases = Dict(
                id => Dict(
                    case_id(n...) => simulate_case(
                        id,
                        n...;
                        thermal_step,
                        burgers_step,
                        burgers_timestep,
                    ) for n in grids_for(id)
                ) for id in keys(expected_cases(family))
            )
            write_fields(
                joinpath(stage, family * ".h5"),
                cases,
                source_metadata(run_id, family),
            )
        end
        read_pair(stage)
    end
    println("N07の両HDF5を保存しました。analyze.jl、plot.jlを再実行してください。")
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
