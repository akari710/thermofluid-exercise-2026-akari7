# N03: 一次元の熱拡散を計算し、境界条件による温度と熱量の違いを調べる。
# 実装する箇所: 内部点更新、断熱端更新、熱量積分のTODO 3か所。
# 読む順: 初期条件 → 境界 → 1ステップ → 熱量 → simulate → main。
# 入力検証・収束診断・描画・保存は provided_support.jl に配布済み。
module N03Diffusion

include("provided_support.jl")

export initial_condition, analytic_solution, apply_boundary!, diffusion_step!, thermal_content,
    simulate, main, stability_experiment

"""
    initial_condition(x; initial=:pulse, boundary=:fixed)

無次元領域 `[0, 2]` の各座標に対応する初期温度を返す。

# 引数

- `x`: 初期温度を与える格子座標の一次元配列。各要素は有限値とする。
- `initial`: 矩形の温度分布を使う `:pulse`、または解析モードを使う `:mode`。
- `boundary`: 固定温度の `:fixed`、または断熱の `:insulated`。

座標が非有限値、または条件の指定が不正な場合は `ArgumentError` を投げる。

# 返り値

`x` の各座標に対応する初期温度の配列。入力の `x` は変更しない。
"""
function initial_condition(x; initial = :pulse, boundary = :fixed)
    # 初期条件と境界の指定、座標の有限性を確認する。
    validate_boundary(boundary)
    validate_initial(initial)
    all(isfinite, x) || throw(ArgumentError("座標を有限値にしてください"))

    # 解析モードは時刻0の解を使い、矩形分布は各座標で温度を決める。
    initial == :mode && return analytic_solution(x, 0.; boundary)
    return [0.5 <= xi <= 1.0 ? 1.0 : 0.0 for xi in x]
end

"""
    analytic_solution(x, t; boundary=:fixed, diffusivity=0.1)

格子収束の確認に使う解析モードの温度を返す。

# 引数

- `x`: 解析解を評価する格子座標の一次元配列。各要素は有限値とする。
- `t`: 有限な非負の時刻。
- `boundary`: 固定温度の `:fixed`、または断熱の `:insulated`。
- `diffusivity`: 有限な正の拡散係数。

値の範囲や境界条件が不正な場合は `ArgumentError` を投げる。

# 返り値

`x` の各座標に対応する時刻 `t` の解析解の配列。入力の `x` は変更しない。
境界条件に適合するモードの解であり、矩形初期値 `initial=:pulse` に対する解ではない。
"""
function analytic_solution(x, t::Real; boundary = :fixed, diffusivity::Real = 0.1)
    # 解析解を評価できる境界条件、拡散係数、時刻、座標を確認する。
    validate_boundary(boundary)
    positive_finite(diffusivity, "diffusivity")
    isfinite(t) && t >= 0 || throw(ArgumentError("tは有限な非負値です"))
    all(isfinite, x) || throw(ArgumentError("座標を有限値にしてください"))

    # 時間による振幅の減衰と、境界に適合する空間分布を組み合わせる。
    amplitude = exp(-diffusivity * (pi / 2)^2 * t)
    return boundary == :fixed ? amplitude .* sin.(pi .* x ./ 2) :
        0.5 .+ 0.5 .* amplitude .* cos.(pi .* x ./ 2)
end

"""
    apply_boundary!(u_new, u_old, r; boundary=:fixed)

新しい温度配列の両端へ、固定温度または断熱の境界条件を適用する。

# 引数

- `u_new`: 両端を書き換える次時刻の温度配列。内部点はこの関数では変更しない。
- `u_old`: 参照する旧時刻の温度配列。全要素を有限値とし、この関数では変更しない。
- `r`: 有限な正の実効Fourier数。
- `boundary`: 両端を温度0にする `:fixed`、または断熱にする `:insulated`。

新旧配列は1始まり・同長・3点以上の浮動小数ベクトルとし、互いに重ならないようにする。
不正な条件は、境界を書き換える前に `ArgumentError` で拒否する。

# 返り値

境界を更新した `u_new` そのもの。新しい配列は作らない。

# 受講生のToDo

断熱端を、旧温度と鏡映ゴースト点を用いて更新する処理を実装する。
固定温度の処理は提供済み。
配布時は、断熱条件を選ぶと `未実装 N03: apply_boundary!` のエラーで停止する。
"""
function apply_boundary!(u_new, u_old, r::Real; boundary = :fixed)
    # 境界を書き込む前に、新旧配列と実効Fourier数を確認する。
    validate_buffers(u_new, u_old, boundary)
    positive_finite(r, "r")

    # 境界条件に応じて両端だけを更新する。
    if boundary == :fixed
        u_new[1] = 0.0
        u_new[end] = 0.0
    else
        # TODO(N03): 旧配列を用い、半セルに対応する係数で両端を更新する。
        error("未実装 N03: apply_boundary!")
    end

    return u_new
end

"""
    diffusion_step!(u_new, u_old, dt, dx, diffusivity; boundary=:fixed)

空間中心差分と陽Euler法で温度を1ステップ進め、両端の境界条件も適用する。

# 引数

- `u_new`: 次時刻の温度を書き込む浮動小数ベクトル。
- `u_old`: 参照する旧時刻の温度ベクトル。全要素を有限値とし、この関数では変更しない。
- `dt`: 有限な正の時間刻み。
- `dx`: 有限な正の格子幅。
- `diffusivity`: 有限な正の拡散係数。
- `boundary`: 固定温度の `:fixed`、または断熱の `:insulated`。

新旧配列は1始まり・同長・3点以上とし、互いに重ならないようにする。
不正な条件や非有限・非正の実効Fourier数は、更新前に `ArgumentError` で拒否する。
実効Fourier数が安定上限を超えていても拒否せず、不安定条件の実験に使える。

# 返り値

実装後は、内部点と境界を更新した `u_new` そのものを返す。新しい配列は作らない。

# 受講生のToDo

旧温度だけを参照して、各内部点の次時刻の温度を求める処理を実装する。
境界処理は、末尾の `apply_boundary!` の呼出しから行う。
配布時は、内部点の更新で `未実装 N03: diffusion_step!` のエラーにより停止する。
"""
function diffusion_step!(
    u_new, u_old, dt::Real, dx::Real, diffusivity::Real; boundary = :fixed,
)
    # 配列の独立性と、刻み・格子幅・拡散係数を確認する。
    validate_buffers(u_new, u_old, boundary)
    for (value, name) in ((dt, "dt"), (dx, "dx"), (diffusivity, "diffusivity"))
        positive_finite(value, name)
    end

    # 更新に使う実効Fourier数を求め、表現可能な正値か確認する。
    r = diffusivity * dt / dx^2
    positive_finite(r, "実効Fourier数")

    # 内部点はすべて同じ旧配列から更新する。
    for i in 2:length(u_old)-1
        # TODO(N03): 旧配列だけから内部点の次時刻値を計算する。
        error("未実装 N03: diffusion_step!")
    end

    # 内部点の後で両端を更新し、次時刻の温度配列を完成させる。
    apply_boundary!(u_new, u_old, r; boundary)
    return u_new
end

"""
    thermal_content(u, dx)

温度配列 `u` を台形則で積分した無次元熱量 `H` を返す。

# 引数

- `u`: 1始まり・3点以上の実数ベクトル。全要素を有限値とし、この関数では変更しない。
- `dx`: 有限な正の格子幅。

不正な条件は `ArgumentError` で拒否する。

# 返り値

実装後は、台形則で求めた熱量を実数で返す。

# 受講生のToDo

端点に内部点の半分の重みを付け、台形則の積分値を返す処理を実装する。
配布時は `未実装 N03: thermal_content` のエラーで停止する。
"""
function thermal_content(u, dx::Real)
    # 積分する温度配列と格子幅を確認する。
    u isa AbstractVector{<:Real} && length(u) >= 3 ||
        throw(ArgumentError("熱量配列は3点以上です"))
    Base.require_one_based_indexing(u)
    all(isfinite, u) || throw(ArgumentError("温度を有限値にしてください"))
    positive_finite(dx, "dx")

    # TODO(N03): 台形則の積分を返す。
    error("未実装 N03: thermal_content")
end

"""
    simulate(; boundary=:fixed, nx=81, diffusivity=0.1, fo=0.4,
               t_final=1.0, initial=:pulse)

無次元領域 `[0, 2]` の熱拡散を指定時刻まで計算し、温度と診断量を返す。

# 引数

- `boundary`: 固定温度の `:fixed` または断熱の `:insulated`。
- `nx`: 両端を含む3点以上の格子点数。`Bool` は不可。
- `diffusivity`: 有限な正の拡散係数。
- `fo`: 要求する有限な正のFourier数。最終時刻に合わせて実効値を調整する。
- `t_final`: 有限な正の終了時刻。
- `initial`: 矩形分布の `:pulse` または格子収束用の `:mode`。

値の範囲や条件の指定が不正な場合は `ArgumentError` を投げる。
安定条件超過や負の温度は実験のため許容する。
計算中に非有限値が生じた場合は、ステップ、時刻、実効Fourier数を示して停止する。
実行には `apply_boundary!`、`diffusion_step!`、`thermal_content` の実装が必要。
未実装の関数へ到達した場合は、その関数の未実装エラーで停止する。

# 返り値

次のフィールドを持つ `NamedTuple`。温度の全時刻の履歴は保持しない。

- `x`, `u0`, `u`: 格子座標、初期温度、最終温度のベクトル。初期温度と最終温度は独立した配列。
- `dx`, `dt`, `steps`: 実際に使った格子幅、時間刻み、ステップ数。
- `requested_fo`, `fo`: 要求したFourier数と、終了時刻に合わせた実効Fourier数。
- `t_final`, `diffusivity`: 終了時刻と拡散係数。
- `times`, `heat_history`: 初期時刻を含む各時刻と、それに対応する熱量のベクトル。
- `initial_heat`, `final_heat`, `heat_change`: 初期熱量、最終熱量、最終値から初期値を引いた変化量。
- `initial_minimum`, `initial_maximum`: 初期温度の最小値と最大値。
- `minimum`, `maximum`: 初期時刻を含む全時刻を通した温度の最小値と最大値。
"""
function simulate(;
    boundary = :fixed, nx::Integer = 81, diffusivity::Real = 0.1,
    fo::Real = 0.4, t_final::Real = 1.0, initial = :pulse,
)
    # 計算条件を確認し、両端を含む格子と初期温度を用意する。
    validate_simulation_inputs(boundary, nx, diffusivity, fo, t_final, initial)
    dx = 2.0 / (nx - 1)
    x = [j * dx for j in 0:nx-1]
    u0 = initial_condition(x; initial, boundary)

    # 正弦の浮動小数丸めも含め、固定端を正確に0へ揃える。
    if boundary == :fixed
        u0[1] = 0.
        u0[end] = 0.
    end

    # 指定Fourier数から刻みを求め、終了時刻へ整数ステップで到達させる。
    nominal_dt = fo * dx^2 / diffusivity
    positive_finite(nominal_dt, "nominal_dt")
    count = t_final / nominal_dt
    isfinite(count) && count < typemax(Int) ||
        throw(ArgumentError("ステップ数が表現範囲を超えています"))
    steps = max(1, ceil(Int, count))
    dt = t_final / steps
    effective_fo = diffusivity * dt / dx^2
    positive_finite(effective_fo, "実効Fourier数")

    # 初期温度を保持し、時間更新用の独立した新旧バッファを用意する。
    u_old = copy(u0)
    u_new = similar(u_old)

    # 初期熱量と極値を記録し、熱量の時刻履歴を用意する。
    initial_heat = thermal_content(u0, dx)
    initial_minimum, initial_maximum = extrema(u0)
    low, high = initial_minimum, initial_maximum
    ensure_finite((initial_heat, low, high), 0, 0., effective_fo)
    times = collect(range(0., t_final; length = steps + 1))
    heat_history = Vector{Float64}(undef, steps + 1)
    heat_history[1] = initial_heat

    # 温度を進め、各時刻の熱量と計算全体の温度範囲を記録する。
    for step in 1:steps
        # 旧温度から次時刻の温度を求め、発散による非有限値を検出する。
        diffusion_step!(u_new, u_old, dt, dx, diffusivity; boundary)
        ensure_finite(u_new, step, times[step+1], effective_fo)

        # 次時刻の熱量と、これまでの温度の極値を記録する。
        heat = thermal_content(u_new, dx)
        ensure_finite((heat,), step, times[step+1], effective_fo)
        heat_history[step+1] = heat
        low = min(low, minimum(u_new))
        high = max(high, maximum(u_new))

        # 計算済みの新温度を次の旧温度にし、空いた配列を再利用する。
        u_old, u_new = u_new, u_old
    end

    # 最終熱量と初期値からの変化をまとめる。
    final_heat = heat_history[end]
    heat_change = final_heat - initial_heat
    ensure_finite((heat_change,), steps, t_final, effective_fo)
    return (;
        x, u0, u = u_old, dx, dt, steps, fo = effective_fo, requested_fo = fo, t_final,
        diffusivity, times, heat_history, initial_heat, final_heat, heat_change,
        initial_minimum, initial_maximum, minimum = low, maximum = high,
    )
end

"""
    main(; output_dir=DEFAULT_OUTPUT_DIR)

両境界の標準計算と解析モードの格子収束を調べ、公式4出力を保存する。

# 引数

- `output_dir`: 出力先のパス。既定値はこの課題の `results/`。

実行には `apply_boundary!`、`diffusion_step!`、`thermal_content` の実装が必要。
未実装の関数へ到達した場合は、その関数の未実装エラーで停止する。
保存処理の失敗時の動作は [`save_staged`](@ref) を参照する。

# 返り値

`(; fixed, insulated, convergence)` という `NamedTuple`。
`fixed` と `insulated` は各境界条件の `simulate` の結果、`convergence` は格子収束の診断辞書。
出力先に温度比較・熱量履歴・収束のPNGと `summary.toml` を保存する。
"""
function main(; output_dir::AbstractString = DEFAULT_OUTPUT_DIR)
    # 固定温度・断熱の標準計算と、解析解に対する格子収束を求める。
    fixed = simulate()
    insulated = simulate(; boundary = :insulated)
    convergence = convergence_results()

    # 診断と図を生成・検査してから、指定された出力先へ保存する。
    write_outputs(output_dir, fixed, insulated, convergence)
    println("N03の出力を書き込みました: $(abspath(output_dir))")
    return (; fixed, insulated, convergence)
end

"""
    stability_experiment(; output_dir, fo=0.6, t_final=0.1)

固定温度で安定・不安定条件を同じ時刻まで進め、比較図と診断を保存する任意実験。

# 引数

- `output_dir`: 出力先のパス。必須で指定し、公式出力と区別できる保存先を選ぶ。
- `fo`: 不安定計算に要求する、有限で0.5より大きいFourier数。
- `t_final`: 両計算で共通の有限な正の終了時刻。

安定計算では要求Fourier数0.4を使う。
終了時刻への調整後も、不安定計算の実効Fourier数が0.5を超える必要がある。
不正な条件は `ArgumentError`、計算中の非有限値はエラーで停止する。
実行にはN03のTODO 3か所の実装が必要。
保存処理の失敗時の動作は [`save_staged`](@ref) を参照する。

# 返り値

`(; stable, unstable)` という `NamedTuple`。両フィールドは `simulate` の結果。
出力先に `stability-comparison.png` と `summary.toml` を保存する。
"""
function stability_experiment(;
    output_dir::AbstractString, fo::Real = 0.6, t_final::Real = 0.1,
)
    # 指定値が不安定条件になることを確認する。
    isfinite(fo) && fo > 0.5 ||
        throw(ArgumentError("比較用foは有限で0.5より大きくしてください"))

    # 安定条件と不安定条件を、同じ初期値・境界・終了時刻で比較する。
    stable = simulate(; fo = 0.4, t_final)
    unstable = simulate(; fo, t_final)
    unstable.fo > 0.5 || throw(ArgumentError(
        "最終時刻への調整後の実効foが0.5以下です。安定条件超過の比較になるよう条件を調整してください",
    ))

    # 任意実験の比較図と診断を保存する。
    write_stability_outputs(output_dir, stable, unstable)
    return (; stable, unstable)
end

# このファイルを直接実行したときだけ、公式出力を作る。
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
