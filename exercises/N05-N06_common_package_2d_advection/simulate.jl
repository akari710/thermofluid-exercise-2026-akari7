# N05/N06: 格子と保存時刻の準備 → 数値更新 → HDF5保存。数値TODOはsrc/を使います。

module N06Simulation
using ThermofluidExercise
include("provided_support.jl")

"""
    simulate_case(
        nx,
        ny;
        times = SAVE_TIMES,
        cx = 1.,
        cy = 0.5,
        safety = 0.8,
        t_final = 1.,
        stepper = ThermofluidExercise.N06.advection_step!,
    )

保存区間ごとに刻みを合わせ、2次元移流の場を記録する。

# 引数

- `nx`: x方向の格子点数。
- `ny`: y方向の格子点数。
- `times`: 先頭0、末尾t_finalの有限な狭義増加の保存時刻列。
- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。
- `safety`: 安定条件に掛ける安全係数。
- `t_final`: 計算の最終時刻。
- `stepper`: 旧場を保持して新場を更新する数値関数。

# 返り値

格子・座標・時刻・u(nx,ny,nt)・区間別刻みとステップ数の辞書。非有限値ではケース・ステップ・時刻・CFLを示して停止する。
"""
function simulate_case(
    nx,
    ny;
    times = SAVE_TIMES,
    cx = 1.,
    cy = 0.5,
    safety = 0.8,
    t_final = 1.,
    stepper = ThermofluidExercise.N06.advection_step!,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        nx isa Integer &&
            ny isa Integer &&
            !(nx isa Bool || ny isa Bool) &&
            nx >= 3 &&
            ny >= 3,
        "nx,nyは3以上の整数です",
    )
    time = schedule(times, t_final)

    # 座標と格子幅を用意し、配列の添字との対応をそろえる。
    dx = 2 / nx
    dy = 1 / ny

    # 安定上限と到達する時刻に合わせて時間刻みを決める。
    dt_max = ThermofluidExercise.N06.stable_timestep(cx, cy, dx, dy; safety)

    # 座標と格子幅を用意し、配列の添字との対応をそろえる。
    x = collect((0:(nx - 1)) .* dx)
    y = collect((0:(ny - 1)) .* dy)

    # 初期値を保持し、更新に使う作業用配列を用意する。
    a = exact_field(x, y, 0.; cx, cy)
    b = similar(a)
    U = Array{Float64}(undef, nx, ny, length(time))
    U[:, :, 1] = a
    dts = Float64[]
    counts = Int[]

    # 旧場から更新し、診断・境界処理を終えてから次の反復へ進む。
    for k in 2:length(time)
        steps, dt = ThermofluidExercise.fit_timestep(dt_max, time[k] - time[k - 1])
        push!(dts, dt)
        push!(counts, steps)
        for step in 1:steps
            stepper(b, a, dt, dx, dy; cx, cy)
            all(isfinite, b) || error(
                "非有限値: case=$(case_id(nx,ny)), step=$step, t=$(time[k - 1] + step * dt), CFL=$(cx * dt / dx + cy * dt / dy)",
            )
            a, b = b, a
        end
        U[:, :, k] = a
    end
    Dict{String,Any}(
        "nx" => nx,
        "ny" => ny,
        "dx" => dx,
        "dy" => dy,
        "dt_max" => dt_max,
        "x" => x,
        "y" => y,
        "time" => time,
        "u" => U,
        "segment_dt" => dts,
        "segment_steps" => counts,
    )
end

"""
    main(;
        output_dir = DEFAULT_OUTPUT_DIR,
        grids = OFFICIAL_GRIDS,
        times = SAVE_TIMES,
        stepper = ThermofluidExercise.N06.advection_step!,
        publish_options...,
    )

公式格子の計算結果を検証してHDF5へ保存する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `grids`: 計算する(nx, ny)の組の列。
- `times`: 先頭0、末尾t_finalの有限な狭義増加の保存時刻列。
- `stepper`: 旧場を保持して新場を更新する数値関数。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

fields.h5のパス。一時生成と読取り検証が成功した場合だけ復元付きで反映する。
"""
function main(;
    output_dir = DEFAULT_OUTPUT_DIR,
    grids = OFFICIAL_GRIDS,
    times = SAVE_TIMES,
    stepper = ThermofluidExercise.N06.advection_step!,
    publish_options...,
)
    # 一時領域で処理を完了してから、公式出力を反映する。
    staged(output_dir, (FIELD_NAME,); publish_options...) do stage
        cases = Dict(
            case_id(nx, ny) => simulate_case(nx, ny; times, stepper) for (nx, ny) in grids
        )
        write_fields(joinpath(stage, FIELD_NAME), cases)
    end
    println("N06 fields.h5を保存しました。analyze.jl、plot.jlを順に再実行してください。")
    joinpath(output_dir, FIELD_NAME)
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
