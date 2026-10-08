# N01: 初期条件、風上差分、中心差分のTODO 3か所を実装する。
# 読む順: 初期条件 → 1ステップ更新 → 境界条件 → simulate → 診断 → main。
# 自分の確認はtests.jl、記録はlearning_log.md。実行結果はresults/へ保存する。
# 入力検証・描画・保存の詳細はprovided_support.jlを参照する。

module N01LinearAdvection

include("provided_support.jl")

export rectangular_initial_condition, apply_boundary!, upwind_step!, centered_step!
export simulate, write_summary, make_plots, main

# === 学生が実装する3つの関数 ===

"""
    rectangular_initial_condition(
        x::AbstractVector{<:Real};
        base::Real = 1.0,
        plateau::Real = 2.0,
        plateau_start::Real = 0.5,
        plateau_end::Real = 1.0,
    )

矩形領域の内側と外側で値を分けた初期分布を作る。

# 引数

- `x`: 空でなく、有限値が狭義単調増加する実数の座標ベクトル。入力は変更しない。
- `base`: 矩形領域の外側の値。
- `plateau`: 矩形領域の内側の値。
- `plateau_start`: 矩形領域の開始座標。
- `plateau_end`: 矩形領域の終了座標。

4パラメータは有限実数で、first(x)からlast(x)の間に開始・終了の順で置く。不正入力はArgumentError。

# 返り値

実装後は `x` と同じ長さの新しい配列。座標は変更しない。

# 受講生のToDo

座標ごとに、両端を含む矩形区間の内側はplateau、外側はbaseとする初期分布を新しい配列に作る。
入力検証は提供済み。F02の配列処理とN01の初期条件を参照する。配布状態では検証後に未実装エラーで停止する。
"""
function rectangular_initial_condition(
    x::AbstractVector{<:Real};
    base::Real = 1.0,
    plateau::Real = 2.0,
    plateau_start::Real = 0.5,
    plateau_end::Real = 1.0,
)
    validate_initial_condition_inputs(x, base, plateau, plateau_start, plateau_end)

    #=
    TODO(N01):
    矩形状の初期分布を実装する。
    各座標xiがplateau_start <= xi <= plateau_endを満たすか判定する。
    =#
    error("未実装 N01: rectangular_initial_condition")
end

"""
    upwind_step!(u_new, u_old, c::Real, dt::Real, dx::Real)

正速度の風上差分と陽Euler法で内部点を更新する。

# 引数

- `u_new`: 内部点を書き込むベクトル。u_oldと記憶領域を共有しないものを用意する。
- `u_old`: 同長・3点以上の有限な旧ベクトル。更新中は変更しない。
- `c`: 有限な正の移流速度。
- `dt`: 有限な正の時間刻み。
- `dx`: x方向の有限な正の格子幅。

提供検証は同じオブジェクトの新旧バッファを拒否する。CFL上限はこの関数では拒否しない。

# 返り値

実装後は更新した `u_new`。両端は保持し、境界条件は呼出し側で適用する。

# 受講生のToDo

提供ループの内部点で、同じ旧ベクトルの現在点と左隣から正速度の風上更新を求める。
入力検証、CFL数courant、旧場のコピーは提供済み。両端は旧値のまま返し、呼出し側が境界を適用する。
N01授業の風上差分と陽Euler法を参照する。配布状態では内部点の未実装エラーで停止する。
"""
function upwind_step!(u_new, u_old, c::Real, dt::Real, dx::Real)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_step_inputs(u_new, u_old, c, dt, dx)
    # courantは1ステップで進む格子幅の割合を表すCFL数です。
    courant = c * dt / dx
    # copyto!で古い値を複製し、内部を更新する間も両端点の値を保ちます。
    copyto!(u_new, u_old)
    for i in 2:(length(u_old) - 1)
        #=
        TODO(N01):
        風上差分と陽Euler法による1ステップの更新を実装する。
        授業ページの式を、iとi - 1の添字を使ってコードへ写す。
        =#
        error("未実装 N01: upwind_step!")
    end
    return u_new
end

"""
    centered_step!(u_new, u_old, c::Real, dt::Real, dx::Real)

中心差分と陽Euler法で内部点を更新し、不安定化を比較する。

# 引数

- `u_new`: 内部点を書き込むベクトル。u_oldと記憶領域を共有しないものを用意する。
- `u_old`: 同長・3点以上の有限な旧ベクトル。更新中は変更しない。
- `c`: 有限な正の移流速度。
- `dt`: 有限な正の時間刻み。
- `dx`: x方向の有限な正の格子幅。

提供検証は同じオブジェクトの新旧バッファを拒否する。CFL上限はこの関数では拒否しない。

# 返り値

実装後は `u_new` 自体。両端は旧値を保ち、境界条件は呼出し側で適用する。`u_old` は保持する。

# 受講生のToDo

提供ループの内部点で、同じ旧ベクトルの左右隣から中心差分と陽Euler法の更新を求める。
入力検証、CFL数courant、旧場のコピーは提供済み。両端は旧値のまま返し、呼出し側が境界を適用する。
N01授業の中心差分と安定性の比較を参照する。配布状態では内部点の未実装エラーで停止する。
"""
function centered_step!(u_new, u_old, c::Real, dt::Real, dx::Real)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_step_inputs(u_new, u_old, c, dt, dx)
    # courantは1ステップで進む格子幅の割合を表すCFL数です。
    courant = c * dt / dx
    # copyto!で古い値を複製し、内部を更新する間も両端点の値を保ちます。
    copyto!(u_new, u_old)
    for i in 2:(length(u_old) - 1)
        #=
        TODO(N01):
        中心差分と陽Euler法による1ステップの更新を実装する。
        授業ページの式を、i - 1とi + 1の添字を使ってコードへ写す。
        =#
        error("未実装 N01: centered_step!")
    end
    return u_new
end

# === 境界条件と時間発展の流れ ===

"""
    apply_boundary!(u::AbstractVector{<:Real}; left_value::Real = 1.0)

左端を固定値、右端を隣接点の値に設定する。

# 引数

- `u`: 境界を書き換える実数ベクトル。2点以上。
- `left_value`: 左端の有限な固定値。

# 返り値

境界を書き換えた入力 `u`。
"""
function apply_boundary!(u::AbstractVector{<:Real}; left_value::Real = 1.0)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_boundary_inputs(u, left_value)
    u[1] = left_value
    u[end] = u[end - 1]
    return u
end

"""
    simulate(;
        scheme,
        nx::Integer = 81,
        c::Real = 1.0,
        cfl::Real = 0.5,
        t_final::Real = 0.5,
    )

初期条件・更新・境界条件を順に適用し、指定時刻まで計算する。

# 引数

- `scheme`: 使用する差分法の識別子。
- `nx`: x方向の格子点数。
- `c`: 正の移流速度。
- `cfl`: 指定するCFL数。
- `t_final`: 計算の最終時刻。

共通の入力検証後、初期条件と選択した差分TODOを呼ぶ。未実装エラーは呼出し元へ伝わる。

# 返り値

`x, u0, u, dx, dt, steps, cfl, minimum, maximum` を持つ `NamedTuple`。初期配列と最終配列は独立する。
"""
function simulate(;
    scheme,
    nx::Integer = 81,
    c::Real = 1.0,
    cfl::Real = 0.5,
    t_final::Real = 0.5,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_simulation_inputs(scheme, nx, c, cfl, t_final)


    # 座標と格子幅を用意し、配列の添字との対応をそろえる。
    x = collect(range(0.0, 2.0; length = nx))
    dx = x[2] - x[1]

    # 安定上限と到達する時刻に合わせて時間刻みを決める。
    nominal_dt = cfl * dx / c
    steps = ceil(Int, t_final / nominal_dt)
    dt = t_final / steps
    actual_cfl = c * dt / dx


    # 初期値を保持し、更新に使う作業用配列を用意する。
    u0 = rectangular_initial_condition(x)
    u_old = copy(u0)
    u_new = similar(u_old)
    # `scheme`に応じて、時間ループで呼ぶ1ステップ関数を選びます。
    step! = scheme === :upwind ? upwind_step! : centered_step!

    # 各ステップで新しい値を計算し、境界条件を適用してから二つのバッファを交換する。
    for _ in 1:steps
        step!(u_new, u_old, c, dt, dx)
        apply_boundary!(u_new)
        u_old, u_new = u_new, u_old
    end


    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return (
        x = x,
        u0 = u0,
        u = u_old,
        dx = dx,
        dt = dt,
        steps = steps,
        cfl = actual_cfl,
        minimum = minimum(u_old),
        maximum = maximum(u_old),
    )
end

"""
    summary_section(scheme::String, result)

初期値の範囲からの超過量を含め、1手法の診断をまとめる。

# 引数

- `scheme`: 使用する差分法の識別子。
- `result`: シミュレーションの結果と診断量。

# 返り値

scheme・CFL・刻み・極値・overshoot/undershootとその判定を持つ辞書。
"""
function summary_section(scheme::String, result)
    initial_minimum, initial_maximum = extrema(result.u0)
    overshoot = max(result.maximum - initial_maximum, 0.0)
    undershoot = max(initial_minimum - result.minimum, 0.0)
    tolerance = 100eps(Float64) * max(abs(initial_minimum), abs(initial_maximum), 1.0)

    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return Dict(
        "scheme" => scheme,
        "cfl" => result.cfl,
        "dt" => result.dt,
        "steps" => result.steps,
        "minimum" => result.minimum,
        "maximum" => result.maximum,
        "overshoot" => overshoot,
        "undershoot" => undershoot,
        "overshoot_occurred" => overshoot > tolerance,
        "undershoot_occurred" => undershoot > tolerance,
    )
end

"""
    main(;
        output_dir::AbstractString = DEFAULT_OUTPUT_DIR,
        nx::Integer = 81,
        c::Real = 1.0,
        cfl::Real = 0.5,
        t_final::Real = 0.5,
    )

風上と中心差分を比較し、診断TOMLと2つの図を保存する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `nx`: x方向の格子点数。
- `c`: 正の移流速度。
- `cfl`: 指定するCFL数。
- `t_final`: 計算の最終時刻。

実行にはN01の3つのTODO実装が必要。計算やファイル保存の失敗は呼出し元へ伝わる。

# 返り値

`upwind, centered, summary_path, plot_paths` を持つ `NamedTuple`。指定ディレクトリへ成果物を書き出す。
"""
function main(;
    output_dir::AbstractString = DEFAULT_OUTPUT_DIR,
    nx::Integer = 81,
    c::Real = 1.0,
    cfl::Real = 0.5,
    t_final::Real = 0.5,
)
    upwind = simulate(; scheme = :upwind, nx, c, cfl, t_final)
    centered = simulate(; scheme = :centered, nx, c, cfl, t_final)
    summary_path = write_summary(output_dir, upwind, centered)
    plot_paths = make_plots(output_dir, upwind, centered)
    println("N01の出力を書き込みました: $(abspath(output_dir))")

    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return (; upwind, centered, summary_path, plot_paths)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
