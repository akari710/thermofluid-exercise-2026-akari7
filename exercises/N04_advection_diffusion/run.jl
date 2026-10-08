# N04: SELECTED_MODELを選び、選択流束・合成刻み・周期更新の3か所を実装する。
# 読む順: 流束 → 刻み → 1ステップ → simulate → main。入力検証・解析解・描画は提供済み。
module N04AdvectionDiffusion
const SELECTED_MODEL = :unselected
include("provided_support.jl")
export SELECTED_MODEL,
    advective_flux,
    stable_timestep,
    advection_diffusion_step!,
    initial_condition,
    analytic_solution,
    conserved_integral,
    simulate,
    main

"""
    advective_flux(u; model = SELECTED_MODEL, speed = 1.0)

選択した線形または非線形モデルの保存形流束を評価する。

# 引数

- `u`: 有限な実数スカラー。非線形モデルでは非負とする。
- `model`: `:linear` または `:nonlinear` の選択モデル。
- `speed`: 線形モデルでは有限な正の実数。非線形モデルでは場の値を速度として使う。

共通検証で不正入力を拒否する。SELECTED_MODELで選択したモデルのTODOを実装する。

# 返り値

実装後はスカラーの流束。

# 受講生のToDo

SELECTED_MODELで選んだ枝だけに、N01の線形流束またはN02のBurgers流束を移す。スカラー1点の入力を扱う。
モデル・有限性・速度の検証は提供済み。配布状態の `:unselected` はモデル検証で拒否され、選択済みの枝も未完成なら既存の未実装エラーで停止する。
"""
function advective_flux(u; model = SELECTED_MODEL, speed = 1.0)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_model(model)
    finite_real(u, "u")
    validate_speed(model, speed)
    if model == :linear
        # TODO(N04): 線形を選んだ場合だけ実装する。
        error("未実装 N04: advective_flux (linear)")
    else
        nonnegative_finite(u, "非線形移流のu")
        # TODO(N04): 非線形を選んだ場合だけ実装する。
        error("未実装 N04: advective_flux (nonlinear)")
    end
end

"""
    stable_timestep(max_speed, dx, diffusivity; safety = 0.8)

移流と拡散の合成安定条件に安全係数を掛けた刻みを求める。

# 引数

- `max_speed`: 有限な非負の最大移流速度。
- `dx`: x方向の有限な正の格子幅。
- `diffusivity`: 有限な非負の拡散係数。
- `safety`: 有限な実数で0より大きく1以下の安全係数。

速度と拡散係数を同時に0にはしない。不正入力は提供検証でArgumentError。刻み選択自体はモデル選択に依存しない。

# 返り値

実装後は有限な正の時間刻み。

# 受講生のToDo

移流と拡散を別々に制限するのでなく、両者の寄与を足した安定上限にsafetyを適用して刻みを返す。
入力検証は提供済み。[N04課題「離散化と時間刻み」](https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/N04.html#discretization)を参照する。配布状態では検証後に未実装エラーで停止する。
"""
function stable_timestep(max_speed, dx, diffusivity; safety = 0.8)
    validate_timestep_inputs(max_speed, dx, diffusivity, safety)
    # TODO(N04): 合成条件から刻みを返す。
    error("未実装 N04: stable_timestep")
end

"""
    advection_diffusion_step!(
        u_new,
        u_old,
        dt,
        dx,
        diffusivity;
        model = SELECTED_MODEL,
        advection = true,
        speed = 1.0,
    )

周期の左右隣接を使い、移流と拡散を同じ旧配列から更新する。

# 引数

- `u_new`: 全点を書き込む1始まりの浮動小数ベクトル。旧場と記憶領域を共有しない。
- `u_old`: 同長・3点以上の有限な旧浮動小数ベクトル。変更しない。
- `dt`: 有限な正の時間刻み。移流・拡散の合成安定条件を満たすこと。
- `dx`: x方向の有限な正の格子幅。
- `diffusivity`: 有限な非負の拡散係数。
- `model`: `:linear` または `:nonlinear` の選択モデル。
- `advection`: 移流を含めるかを表す `Bool`。
- `speed`: 線形モデルの移流速度。非線形モデルでは場の値を速度として使う。

モデルと速度・配列・合成安定条件の違反は書込み前に提供検証でArgumentError。配布状態のSELECTED_MODEL=:unselectedも拒否する。

# 返り値

実装後は全点を書き換えた `u_new`。旧配列を保持する。

# 受講生のToDo

提供ループの全点で、周期の左右隣接を使い、移流と拡散を同じ旧ベクトルから評価して新ベクトルへ書く。
新旧は同じ記憶領域を共有せず、更新中の新値を旧値の代わりに読まない。入力・合成安定条件の検証は提供済み。
移流を含める場合は選択モデルのadvective_fluxも必要。advection=falseでは拡散だけを扱う。[N04課題「離散化と時間刻み」](https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/N04.html#discretization)を参照する。配布状態では検証後に未実装エラーで停止する。
"""
function advection_diffusion_step!(
    u_new,
    u_old,
    dt,
    dx,
    diffusivity;
    model = SELECTED_MODEL,
    advection = true,
    speed = 1.0,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_step(u_new, u_old, dt, dx, diffusivity, model, advection, speed)
    for i in eachindex(u_old)
        # TODO(N04): 周期の左右隣接を使い、移流と拡散を同じ旧配列から足す。
        error("未実装 N04: advection_diffusion_step!")
    end
    return u_new
end

"""
    simulate(;
        model = SELECTED_MODEL,
        advection = true,
        nx = 80,
        diffusivity = 0.1,
        speed = 1.0,
        safety = 0.8,
        t_final = 1.0,
        initial = :pulse,
        dt = nothing,
    )

合成条件の刻みを最終時刻に合わせ、周期の移流拡散を計算する。

# 引数

- `model`: `:linear` または `:nonlinear` の選択モデル。
- `advection`: 移流を含めるかを表す `Bool`。
- `nx`: x方向の格子点数。
- `diffusivity`: 拡散係数。
- `speed`: 線形モデルの移流速度。非線形モデルでは場の値を速度として使う。
- `safety`: 安定条件に掛ける安全係数。
- `t_final`: 計算の最終時刻。
- `initial`: 初期条件の識別子。
- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。

指定 `dt` は上限として最終時刻に合わせる。不正入力・合成条件違反は拒否し、非有限値では診断を示して停止する。TODO実装が必要。

# 返り値

`x, u0, u`、モデル・条件、格子/刻み、時刻と積分履歴、初期/最終積分・極値・最大CFL・Fo・合成安定数を持つ `NamedTuple`。
"""
function simulate(;
    model = SELECTED_MODEL,
    advection = true,
    nx = 80,
    diffusivity = 0.1,
    speed = 1.0,
    safety = 0.8,
    t_final = 1.0,
    initial = :pulse,
    dt = nothing,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_simulation(model, advection, nx, diffusivity, speed, safety, t_final, initial)

    # 座標と格子幅を用意し、配列の添字との対応をそろえる。
    dx = 2.0 / nx
    x = [j * dx for j in 0:(nx - 1)]

    # 初期値を保持し、更新に使う作業用配列を用意する。
    u0 = initial_condition(x; model, initial, diffusivity, speed)
    a = maximum_speed(u0, model, advection, speed)

    # 安定上限と到達する時刻に合わせて時間刻みを決める。
    nominal_dt = stable_timestep(a, dx, diffusivity; safety)
    positive_finite(nominal_dt, "nominal_dt")
    if !isnothing(dt)
        positive_finite(dt, "dt")
        stability_numbers(a, dt, dx, diffusivity)
        nominal_dt = dt
    end
    count = t_final / nominal_dt
    isfinite(count) && count < typemax(Int) - 1 ||
        throw(ArgumentError("ステップ数が表現範囲を超えています"))
    steps = max(1, ceil(Int, count))
    dt = t_final / steps
    positive_finite(dt, "実効dt")
    cfl, fo = stability_numbers(a, dt, dx, diffusivity)

    # 初期値を保持し、更新に使う作業用配列を用意する。
    u_old = copy(u0)
    u_new = similar(u_old)
    times = collect(range(0., t_final; length = steps + 1))
    integral_history = Vector{Float64}(undef, steps + 1)

    # 初期値を基準にして、保存量や停止条件の診断を準備する。
    initial_integral = conserved_integral(u0, dx)
    integral_history[1] = initial_integral
    low, high = extrema(u0)
    max_cfl = cfl

    # 旧場から更新し、診断・境界処理を終えてから次の反復へ進む。
    for step in 1:steps
        a = maximum_speed(u_old, model, advection, speed)
        cfl, fo = stability_numbers(a, dt, dx, diffusivity)
        advection_diffusion_step!(
            u_new,
            u_old,
            dt,
            dx,
            diffusivity;
            model,
            advection,
            speed,
        )
        ensure_finite(u_new, model, step, times[step + 1], cfl, fo)
        integral = dx * sum(u_new)
        ensure_finite((integral,), model, step, times[step + 1], cfl, fo)
        integral_history[step + 1] = integral
        low = min(low, minimum(u_new))
        high = max(high, maximum(u_new))
        a = maximum_speed(u_new, model, advection, speed)
        cfl, fo = stability_numbers(a, dt, dx, diffusivity)
        max_cfl = max(max_cfl, cfl)
        u_old, u_new = u_new, u_old
    end

    # 初期値を基準にして、保存量や停止条件の診断を準備する。
    final_integral = integral_history[end]
    integral_change = final_integral - initial_integral
    max_stability_number = max_cfl + 2fo
    ensure_finite(
        (integral_change, low, high, max_stability_number),
        model,
        steps,
        t_final,
        max_cfl,
        fo,
    )

    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return (;
        x,
        u0,
        u = u_old,
        model,
        advection,
        initial,
        speed,
        diffusivity,
        dx,
        dt,
        steps,
        t_final,
        requested_safety = safety,
        times,
        integral_history,
        initial_integral,
        final_integral,
        integral_change,
        minimum = low,
        maximum = high,
        max_cfl,
        fo,
        max_stability_number,
    )
end

"""
    main(; model = SELECTED_MODEL, output_dir = DEFAULT_OUTPUT_DIR)

移流のみ・拡散のみ・合成の比較と解析モードの格子収束を保存する。

# 引数

- `model`: `:linear` または `:nonlinear` の選択モデル。
- `output_dir`: 公式成果物を書き出すディレクトリ。

# 返り値

`advection_only, diffusion_only, combined, convergence` を持つ `NamedTuple`。`output_dir/<model>/` に図3枚と診断TOMLを保存する。
"""
function main(; model = SELECTED_MODEL, output_dir = DEFAULT_OUTPUT_DIR)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_model(model)
    combined = simulate(; model)
    advection_only = simulate(; model, diffusivity = 0., dt = combined.dt)
    diffusion_only = simulate(; model, advection = false, dt = combined.dt)
    convergence = convergence_results(model)
    write_outputs(output_dir, advection_only, diffusion_only, combined, convergence)
    println("N04の出力を書き込みました: $(abspath(joinpath(output_dir,string(model))))")

    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return (; advection_only, diffusion_only, combined, convergence)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
end
