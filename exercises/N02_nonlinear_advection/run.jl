# N02: 流束、周期左隣添字、流束差分のTODO 3か所を実装する。
# 読む順: 流束 → 添字 → 1ステップ → 境界 → simulate → main。
module N02NonlinearAdvection
include("provided_support.jl")
export burgers_flux,
    periodic_left_index, nonlinear_upwind_step!, apply_boundary!, simulate, main

"""
    burgers_flux(u::Real)

有限なスカラーの保存形Burgers流束を評価する。

# 引数

- `u`: 1点の場を表す有限な実数スカラー。負値も受け付ける。

`u` は有限な実数。条件違反は `ArgumentError`。

# 返り値

実装後はスカラーの流束。流束自体は負の入力にも定義する。

# 受講生のToDo

N02授業の保存形Burgers流束を、配列全体ではなく1点のスカラー値から評価する。負値にも流束を定義する。
有限性検証は提供済み。配布状態では検証後に未実装エラーで停止する。
"""
function burgers_flux(u::Real)
    isfinite(u) || throw(ArgumentError("流束の入力を有限値にしてください"))
    # TODO(N02): 保存形Burgers方程式の流束を返す。
    error("未実装 N02: burgers_flux")
end

"""
    periodic_left_index(i::Integer, n::Integer)

先頭の左隣を末尾にする周期添字を求める。

# 引数

- `i`: 1始まりの現在の添字。
- `n`: 配列の点数。

`i` と `n` はBoolを除く整数で、`n >= 1`、`1 <= i <= n`。不正入力は `ArgumentError`。

# 返り値

実装後は1から `n` までの整数添字。

# 受講生のToDo

先頭の左隣は末尾、それ以外は直前の点とする整数添字を返す。1点の配列でも有効な添字を返す。
入力検証は提供済み。N02の周期境界を参照する。配布状態では検証後に未実装エラーで停止する。
"""
function periodic_left_index(i::Integer, n::Integer)
    (i isa Bool || n isa Bool || n < 1 || !(1 <= i <= n)) &&
        throw(ArgumentError("添字は1 <= i <= n、点数はn >= 1の整数です（Bool不可）"))
    # TODO(N02): 先頭の左隣が末尾になるようにする。
    error("未実装 N02: periodic_left_index")
end

"""
    nonlinear_upwind_step!(u_new, u_old, dt::Real, dx::Real; boundary = :fixed)

旧配列の流束後退差分を使い、陽Euler法で新配列を更新する。

# 引数

- `u_new`: 更新結果を書き込む1始まりの浮動小数ベクトル。旧場と記憶領域を共有しない。
- `u_old`: 同長・3点以上の有限な旧浮動小数ベクトル。変更しない。
- `dt`: 有限な正の時間刻み。CFL超過でも拒否しない。
- `dx`: x方向の有限な正の格子幅。
- `boundary`: `:fixed`（両端の旧値を保つ）または `:periodic`（全点を更新する）。

条件違反は書込み前にArgumentError。更新後の非有限値は提供の検査でエラーになる。

# 返り値

実装後は `u_new`。固定境界では内部点、周期境界では全点を更新し、`u_old` は変更しない。

# 受講生のToDo

提供ループの各点で、旧場の現在点と左隣のBurgers流束の差から更新する。
検証、旧場のコピー、境界別の更新範囲と左隣添字は提供済み。固定境界の両端は旧値を保ち、周期境界は全点を更新する。
新旧は同じ記憶領域を共有せず、全点で更新前の値を読む。N02授業の保存形更新を参照する。配布状態では更新範囲内の未実装エラーで停止する。
"""
function nonlinear_upwind_step!(u_new, u_old, dt::Real, dx::Real; boundary = :fixed)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_step_inputs(u_new, u_old, dt, dx, boundary)
    copyto!(u_new, u_old)
    indices = boundary == :periodic ? (1:length(u_old)) : (2:(length(u_old) - 1))
    for i in indices
        left = boundary == :periodic ? periodic_left_index(i, length(u_old)) : i - 1
        # TODO(N02): 古い値の流束の差で更新する。u_oldは変更しない。
        error("未実装 N02: nonlinear_upwind_step!")
    end
    all(isfinite, u_new) || error("更新後に非有限値があります。更新式を確認してください")
    return u_new
end

"""
    apply_boundary!(u; boundary = :fixed)

固定境界の端点を設定し、周期境界では配列を保持する。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `boundary`: 境界条件の識別子。

配列は3点以上で、境界は `:fixed` または `:periodic`。不正入力は `ArgumentError`。

# 返り値

入力 `u`。固定の場合は左端を1、右端を左隣の値に書き換える。
"""
function apply_boundary!(u; boundary = :fixed)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_boundary(boundary)
    length(u) >= 3 || throw(ArgumentError("境界配列には3点以上が必要です"))
    if boundary == :fixed
        u[1] = 1.0
        u[end] = u[end - 1]
    end
    return u
end

"""
    simulate(;
        boundary = :fixed,
        nx::Integer = 81,
        cfl::Real = 0.5,
        t_final::Real = 1.0,
    )

初期最大速度で刻みを決め、指定時刻までBurgers移流を計算する。

# 引数

- `boundary`: 境界条件の識別子。
- `nx`: x方向の格子点数。
- `cfl`: 指定するCFL数。
- `t_final`: 計算の最終時刻。

固定は81点、周期の公式条件は80点。CFL超過でも計算を続けるが、非有限値になれば停止する。3つのTODO実装が必要。

# 返り値

座標 `x`、独立した初期値 `u0` と最終値 `u`、刻み・ステップ数・最大CFL、極値・超過量・和の変化を持つ `NamedTuple`。
"""
function simulate(;
    boundary = :fixed,
    nx::Integer = 81,
    cfl::Real = 0.5,
    t_final::Real = 1.0,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_simulation_inputs(boundary, nx, cfl, t_final)

    # 座標と格子幅を用意し、配列の添字との対応をそろえる。
    dx = 2.0 / (boundary == :fixed ? nx - 1 : nx)
    x = [j * dx for j in 0:(nx - 1)]
    u0 = [0.5 <= xi <= 1.0 ? 2.0 : 1.0 for xi in x]

    # 安定上限と到達する時刻に合わせて時間刻みを決める。
    nominal_dt = cfl * dx / maximum(abs, u0)
    steps = ceil(Int, t_final / nominal_dt)
    dt = t_final / steps

    # 初期値を保持し、更新に使う作業用配列を用意する。
    u_old = copy(u0)
    u_new = similar(u_old)
    max_cfl = maximum(abs, u_old) * dt / dx

    # 旧場から更新し、診断・境界処理を終えてから次の反復へ進む。
    for _ in 1:steps
        nonlinear_upwind_step!(u_new, u_old, dt, dx; boundary)
        apply_boundary!(u_new; boundary)
        all(isfinite, u_new) || error("更新後に非有限値があります")
        max_cfl = max(max_cfl, maximum(abs, u_new) * dt / dx)
        u_old, u_new = u_new, u_old
    end
    initial_minimum, initial_maximum = extrema(u0)
    low, high = extrema(u_old)

    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return (;
        x,
        u0,
        u = u_old,
        dx,
        dt,
        steps,
        max_cfl,
        t_final,
        initial_minimum,
        initial_maximum,
        minimum = low,
        maximum = high,
        overshoot = max(high - initial_maximum, 0.0),
        undershoot = max(initial_minimum - low, 0.0),
        initial_sum = sum(u0),
        final_sum = sum(u_old),
        sum_change = sum(u_old) - sum(u0),
    )
end

"""
    main(;
        output_dir::AbstractString = DEFAULT_OUTPUT_DIR,
        cfl::Real = 0.5,
        t_final::Real = 1.0,
    )

固定81点・周期80点の公式計算と3成果物を生成する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `cfl`: 指定するCFL数。
- `t_final`: 計算の最終時刻。

# 返り値

`fixed` と `periodic` の計算結果を持つ `NamedTuple`。図2枚と `summary.toml` を保存する。
"""
function main(;
    output_dir::AbstractString = DEFAULT_OUTPUT_DIR,
    cfl::Real = 0.5,
    t_final::Real = 1.0,
)
    fixed = simulate(; boundary = :fixed, nx = 81, cfl, t_final)
    periodic = simulate(; boundary = :periodic, nx = 80, cfl, t_final)
    write_outputs(output_dir, fixed, periodic, cfl, t_final)
    println("N02の出力を書き込みました: $(abspath(output_dir))")

    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return (; fixed, periodic)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
end
