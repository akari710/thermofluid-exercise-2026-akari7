# 入力検証・解析解・診断・描画・保存の提供部分。学生の編集対象ではありません。
using Plots, TOML
const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "results")

"""
    finite_real(v, name)

有限な実数か確認する。

# 引数

- `v`: 検証する値。
- `name`: 検証エラーに表示する項目名。

`x` は有限な実数。不正入力は `ArgumentError`。

# 返り値

検証成功時は `nothing`。 入力は変更しない。
"""
finite_real(v, name) =
    v isa Real && isfinite(v) ? nothing : throw(ArgumentError("$name は有限な実数です"))

"""
    nonnegative_finite(v, name)

有限な非負実数か確認する。

# 引数

- `v`: 検証する値。
- `name`: 検証エラーに表示する項目名。

`x` は有限な非負実数。不正入力は `ArgumentError`。

# 返り値

検証成功時は `true`。 入力は変更しない。
"""
function nonnegative_finite(v, name)
    finite_real(v, name)
    v >= 0 || throw(ArgumentError("$name は非負です"))
end

"""
    positive_finite(v, name)

有限な正の実数か確認する。

# 引数

- `v`: 検証する値。
- `name`: 検証エラーに表示する項目名。

`x` は有限な正の実数。不正入力は `ArgumentError`。

# 返り値

検証成功時は `true`。 入力は変更しない。
"""
function positive_finite(v, name)
    finite_real(v, name)
    v > 0 || throw(ArgumentError("$name は正です"))
end

"""
    validate_model(m)

選択モデルが線形または非線形か確認する。

# 引数

- `m`: 検証するモデル名。

`m` は `:linear` または `:nonlinear`。未選択を含む不正入力は `ArgumentError`。

# 返り値

検証成功時は `nothing`。 入力は変更しない。
"""
validate_model(m) =
    m in (:linear, :nonlinear) ? nothing :
    throw(
        ArgumentError(
            "SELECTED_MODELを:linearまたは:nonlinearに設定してください（model=$m）",
        ),
    )

"""
    validate_speed(model, speed)

線形モデルでは正速度、非線形では有限速度か確認する。

# 引数

- `model`: `:linear` または `:nonlinear` の選択モデル。
- `speed`: 線形モデルの移流速度。非線形モデルでは場の値を速度として使う。

# 返り値

検証成功時は 線形は `true`、非線形は `false`。 入力は変更しない。
"""
function validate_speed(model, speed)
    finite_real(speed, "speed")
    model == :linear && positive_finite(speed, "線形speed")
end

"""
    validate_safety(safety)

安全係数が0より大きく1以下か確認する。

# 引数

- `safety`: 安定条件に掛ける安全係数。

`s` は有限値で `0 < s <= 1`。不正入力は `ArgumentError`。

# 返り値

検証成功時は `true`。 入力は変更しない。
"""
function validate_safety(safety)
    positive_finite(safety, "safety")
    safety <= 1 || throw(ArgumentError("safetyは0より大きく1以下です"))
end

"""
    validate_timestep_inputs(a, dx, d, safety)

速度・格子幅・拡散係数・安全係数を検証する。

# 引数

- `a`: 有限な非負の最大移流速度。
- `dx`: x方向の有限な正の格子幅。
- `d`: 拡散係数。
- `safety`: 安定条件に掛ける安全係数。

速度と拡散係数は有限非負で同時に0にしない。格子幅は有限正値、`0 < safety <= 1`。
不正入力は `ArgumentError`。

# 返り値

検証成功時は `true`。 入力は変更しない。
"""
function validate_timestep_inputs(a, dx, d, safety)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    nonnegative_finite(a, "max_speed")
    positive_finite(dx, "dx")
    nonnegative_finite(d, "diffusivity")
    validate_safety(safety)
    a > 0 ||
        d > 0 ||
        throw(ArgumentError("移流速度とdiffusivityがともに0では刻みを定義できません"))
end

"""
    checked_timestep(dt)

刻みの正値と表現範囲を検証する。

# 引数

- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。

# 返り値

入力 `dt`。
"""
function checked_timestep(dt)
    positive_finite(dt, "合成条件のdt（係数・dxの表現範囲を確認）")
    return dt
end

"""
    maximum_speed(u, model, advection, speed)

移流の有無と選択モデルに応じた最大速度を評価する。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `model`: `:linear` または `:nonlinear` の選択モデル。
- `advection`: 移流を含めるかを表す `Bool`。
- `speed`: 線形モデルの移流速度。非線形モデルでは場の値を速度として使う。

線形モデルは指定速度、非線形モデルは有限な非負の場の値を使う。負の非線形速度は `ArgumentError`。

# 返り値

移流なしは0、線形は `speed`、非線形は場の最大値。非線形の負速度は拒否する。
"""
function maximum_speed(u, model, advection, speed)
    !advection && return 0.
    if model == :nonlinear
        minimum(u) >= 0 ||
            throw(ArgumentError("非線形移流は非負速度の範囲だけに対応します"))
        return maximum(u)
    end
    return speed
end

"""
    stability_numbers(a, dt, dx, d)

CFLとFoを評価して合成安定条件を検証する。

# 引数

- `a`: 有限な非負の最大移流速度。
- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。
- `dx`: x方向の有限な正の格子幅。
- `d`: 拡散係数。

CFLとFoが有限で、合成安定条件 `CFL + 2Fo <= 1` を満たすこと。条件違反は `ArgumentError`。

# 返り値

`(cfl, fo)`。
"""
function stability_numbers(a, dt, dx, d)
    cfl = a * (dt / dx)
    fo = d * (dt / dx) / dx
    all(isfinite, (cfl, fo)) || throw(ArgumentError("CFL・Foが有限ではありません"))
    cfl + 2fo <= 1 + 32eps(Float64) ||
        throw(ArgumentError("合成安定条件CFL+2Fo<=1を超えています: CFL=$cfl, Fo=$fo"))
    return cfl, fo
end

"""
    validate_step(a, b, dt, dx, d, model, advection, speed)

旧配列の保持と合成安定条件に必要な入力を検証する。

# 引数

- `a`: 更新結果を書き込む新ベクトル。旧ベクトルと重ならない。
- `b`: 読み取り専用の旧ベクトル。
- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。
- `dx`: x方向の有限な正の格子幅。
- `d`: 拡散係数。
- `model`: `:linear` または `:nonlinear` の選択モデル。
- `advection`: 移流を含めるかを表す `Bool`。
- `speed`: 線形モデルの移流速度。非線形モデルでは場の値を速度として使う。

新旧は1始まり・同長・3点以上の浮動小数ベクトルで、互いに重ならない。旧場は有限値。
刻みと格子幅は有限正値、拡散係数は有限非負、モデルと速度および合成安定条件も検証する。
不正入力は書込み前に `ArgumentError`。

# 返り値

検証成功時は CFLとFoの2要素タプル。 入力は変更しない。
"""
function validate_step(a, b, dt, dx, d, model, advection, speed)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_model(model)
    validate_speed(model, speed)
    advection isa Bool || throw(ArgumentError("advectionはBoolです"))
    all(u->u isa AbstractVector{<:AbstractFloat}, (a, b)) ||
        throw(ArgumentError("新旧配列は浮動小数の一次元配列です"))
    Base.require_one_based_indexing(a, b)
    length(a) == length(b) >= 3 || throw(ArgumentError("新旧配列は同長で3点以上です"))
    Base.mightalias(a, b) && throw(ArgumentError("新旧配列が重なっています"))
    all(isfinite, b) || throw(ArgumentError("旧配列は有限値にしてください"))
    positive_finite(dt, "dt")
    positive_finite(dx, "dx")
    nonnegative_finite(d, "diffusivity")
    stability_numbers(maximum_speed(b, model, advection, speed), dt, dx, d)
end

"""
    validate_simulation(model, advection, nx, d, speed, safety, t, initial)

計算モデルと格子・係数・時刻・初期条件を検証する。

# 引数

- `model`: `:linear` または `:nonlinear` の選択モデル。
- `advection`: 移流を含めるかを表す `Bool`。
- `nx`: x方向の格子点数。
- `d`: 拡散係数。
- `speed`: 線形モデルの移流速度。非線形モデルでは場の値を速度として使う。
- `safety`: 安定条件に掛ける安全係数。
- `t`: 評価する時刻。
- `initial`: 初期条件の識別子。

モデルを選択し、`advection` はBool、`nx` はBoolを除く3以上の整数、`initial` は `:pulse` または `:mode`。
係数・刻みの安全係数・最終時刻も検証し、不正入力は `ArgumentError`。

# 返り値

検証成功時は `true`。 入力は変更しない。
"""
function validate_simulation(model, advection, nx, d, speed, safety, t, initial)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_model(model)
    validate_speed(model, speed)
    validate_safety(safety)
    advection isa Bool || throw(ArgumentError("advectionはBoolです"))
    nx isa Integer && !(nx isa Bool) && nx >= 3 ||
        throw(ArgumentError("nxは3以上の整数です（Bool不可）"))
    nonnegative_finite(d, "diffusivity")
    positive_finite(t, "t_final")
    initial in (:pulse, :mode) || throw(ArgumentError("initialは:pulseまたは:modeです"))
end

"""
    ensure_finite(values, model, step, t, cfl, fo)

非有限値を検出し、モデル・時刻・安定数を示して停止する。

# 引数

- `values`: 有限性を検査する値の列。
- `model`: `:linear` または `:nonlinear` の選択モデル。
- `step`: 現在のステップ番号。
- `t`: 評価する時刻。
- `cfl`: 指定するCFL数。
- `fo`: 拡散数。

# 返り値

検証成功時は `true`。 入力は変更しない。
"""
function ensure_finite(values, model, step, t, cfl, fo)
    all(isfinite, values) ||
        error("非有限値のため停止: model=$model, step=$step, t=$t, CFL=$cfl, Fo=$fo")
end

"""
    initial_condition(
        x;
        model = SELECTED_MODEL,
        initial = :pulse,
        diffusivity = 0.1,
        speed = 1.0,
    )

周期領域の矩形または解析モードの初期配列を作る。

# 引数

- `x`: x方向の座標。入力は変更しない。
- `model`: `:linear` または `:nonlinear` の選択モデル。
- `initial`: 初期条件の識別子。
- `diffusivity`: 拡散係数。
- `speed`: 線形モデルの移流速度。非線形モデルでは場の値を速度として使う。

座標は有限値、`initial` は `:pulse` または `:mode`。係数とモデルも検証し、不正入力は `ArgumentError`。

# 返り値

`x` に対応する新しい初期配列。
"""
function initial_condition(
    x;
    model = SELECTED_MODEL,
    initial = :pulse,
    diffusivity = 0.1,
    speed = 1.0,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_model(model)
    validate_speed(model, speed)
    nonnegative_finite(diffusivity, "diffusivity")
    all(isfinite, x) || throw(ArgumentError("xは有限値です"))
    initial in (:pulse, :mode) || throw(ArgumentError("initialは:pulseまたは:modeです"))
    initial == :mode && return analytic_solution(x, 0.; model, diffusivity, speed)
    return [0.5 <= xi <= 1. ? 2. : 1. for xi in x]
end

"""
    analytic_solution(x, t; model = SELECTED_MODEL, diffusivity = 0.1, speed = 1.0)

滑らかな解析モードの時刻tの場を評価する。

# 引数

- `x`: x方向の座標。入力は変更しない。
- `t`: 評価する時刻。
- `model`: `:linear` または `:nonlinear` の選択モデル。
- `diffusivity`: 拡散係数。
- `speed`: 線形モデルの移流速度。非線形モデルでは場の値を速度として使う。

入力条件を満たさない場合は `ArgumentError`。

# 返り値

`x` に対応する新しい配列。矩形初期条件の解析解ではない。非線形は平均速度1で、正の拡散係数と非負速度を満たす範囲を求める。
"""
function analytic_solution(x, t; model = SELECTED_MODEL, diffusivity = 0.1, speed = 1.0)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_model(model)
    validate_speed(model, speed)
    nonnegative_finite(diffusivity, "diffusivity")
    nonnegative_finite(t, "t")
    all(isfinite, x) || throw(ArgumentError("xは有限値です"))
    if model == :linear
        return 1 .+ 0.5exp(-diffusivity * pi ^ 2 * t) .* sin.(pi .* (x .- speed * t))
    end
    positive_finite(diffusivity, "非線形解析モードのdiffusivity")
    # max |u-1| at t=0 is 2D*k*A/sqrt(1-A^2), A=0.5.
    diffusivity * pi / sqrt(0.75) <= 1 ||
        throw(ArgumentError("非線形解析モードで負速度を生むdiffusivityです"))
    a = 0.5exp(-diffusivity * pi ^ 2 * t)
    theta = pi .* (x .- t)
    return 1 .+ 2diffusivity * pi * a .* sin.(theta) ./ (1 .+ a .* cos.(theta))
end

"""
    conserved_integral(u, dx)

終点を重複させない周期セルの離散積分を求める。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `dx`: x方向の有限な正の格子幅。

入力条件を満たさない場合は `ArgumentError`。

# 返り値

有限な `dx * sum(u)`。
"""
function conserved_integral(u, dx)
    u isa AbstractVector{<:Real} && length(u) >= 3 ||
        throw(ArgumentError("積分配列は3点以上の実数ベクトルです"))
    Base.require_one_based_indexing(u)
    all(isfinite, u) || throw(ArgumentError("積分配列は有限値です"))
    positive_finite(dx, "dx")
    result = dx * sum(u)

    # 計算や書込みの前に、入力条件をまとめて確認する。
    finite_real(result, "離散積分")
    return result
end

"""
    convergence_results(model)

40・80・160点で滑らかなモードの格子収束を評価する。

# 引数

- `model`: `:linear` または `:nonlinear` の選択モデル。

# 返り値

初期条件名・格子点数・幅・最大誤差・観測次数の辞書。
"""
function convergence_results(model)
    nx = [40, 80, 160]
    runs = [simulate(; model, nx = n, initial = :mode) for n in nx]
    errors = [maximum(abs.(r.u - analytic_solution(r.x, r.t_final; model))) for r in runs]
    for (r, error) in zip(runs, errors)
        ensure_finite((error,), model, r.steps, r.t_final, r.max_cfl, r.fo)
    end
    orders = log2.(errors[1:2] ./ errors[2:3])
    r = last(runs)
    ensure_finite(orders, model, r.steps, r.t_final, r.max_cfl, r.fo)

    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return Dict(
        "initial" => "mode",
        "nx" => nx,
        "dx" => [r.dx for r in runs],
        "errors" => errors,
        "orders" => orders,
    )
end

"""
    summary_section(r)

1条件の有限な診断量と速度の意味をまとめる。

# 引数

- `r`: シミュレーションの結果と診断量。

# 返り値

モデル・係数・刻み・積分・安定数の辞書。
"""
function summary_section(r)
    # 計算値と条件を、保存する診断情報にまとめる。
    section = Dict{String,Any}(
        "nx" => length(r.x),
        "model" => string(r.model),
        "advection" => r.advection,
        "initial" => string(r.initial),
        "speed_role" =>
            r.model == :linear ? "constant advection speed" :
            "unused; velocity is u, analytic mean is 1",
    )
    for key in (
        :speed,
        :diffusivity,
        :dx,
        :dt,
        :steps,
        :t_final,
        :requested_safety,
        :initial_integral,
        :final_integral,
        :integral_change,
        :minimum,
        :maximum,
        :max_cfl,
        :fo,
        :max_stability_number,
    )
        value = getproperty(r, key)
        ensure_finite((value,), r.model, r.steps, r.t_final, r.max_cfl, r.fo)
        section[string(key)] = value
    end
    return section
end

"""
    periodic_display(r, u)

描画専用に周期終点の座標と先頭値を追加する。

# 引数

- `r`: シミュレーションの結果と診断量。
- `u`: 場の値。配列の添字は座標の順に対応する。

# 返り値

拡張した座標と値のタプル。入力配列を変更しない。
"""
periodic_display(r, u) = (vcat(r.x, 2.), vcat(u, u[1]))
const CONDITION_COLORS = ("#0072B2", "#D55E00", "#009E73")

"""
    make_plots(directory, a, b, c, convergence)

3条件の分布・保存誤差・格子収束の公式図を保存する。

# 引数

- `directory`: 出力先ディレクトリ。
- `a`: 移流のみの計算結果。
- `b`: 拡散のみの計算結果。
- `c`: 移流と拡散を含む計算結果。
- `convergence`: 滑らかな解析モードの格子収束結果。

# 返り値

最後の `savefig` の結果。指定先に `comparison.png, conservation.png, convergence.png` を保存する。
"""
function make_plots(directory, a, b, c, convergence)
    p = plot(
        periodic_display(c, c.u0)...;
        label = "Initial (t = 0)",
        color = :gray,
        linestyle = :dash,
        linewidth = 2,
        xlabel = "x (dimensionless)",
        ylabel = c.model == :linear ? "Temperature u (dimensionless)" :
                 "Velocity u (dimensionless)",
        title = "$(c.model), t = $(c.t_final)",
        size = (800, 500),
        ylims = (0.95, 2.05),
        legend = :outerright,
    )
    for ((r, label), color) in zip(
        ((a, "Advection only"), (b, "Diffusion only"), (c, "Combined")),
        CONDITION_COLORS,
    )
        plot!(
            p,
            periodic_display(r, r.u)...;
            label,
            color,
            linestyle = :solid,
            linewidth = 2,
        )
    end
    savefig(p, joinpath(directory, "comparison.png"))
    conservation_title =
        c.model == :linear ? "Temperature integral conservation error" :
        "Velocity integral conservation error"
    p = plot(;
        xlabel = "t (dimensionless)",
        ylabel = "(I - I0) / 1e-14",
        size = (800, 500),
        title = "$conservation_title\n$(c.model): I0 = $(round(c.initial_integral;digits = 5)) (dimensionless)",
        titlefontsize = 12,
        legend = :outerright,
    )
    for ((r, label), color) in zip(
        ((a, "Advection only"), (b, "Diffusion only"), (c, "Combined")),
        CONDITION_COLORS,
    )
        plot!(
            p,
            r.times,
            (r.integral_history .- r.initial_integral) ./ 1e-14;
            label,
            color,
            linestyle = :solid,
            linewidth = 2,
        )
    end
    savefig(p, joinpath(directory, "conservation.png"))

    # 座標と格子幅を用意し、配列の添字との対応をそろえる。
    dx = convergence["dx"]
    errors = convergence["errors"]
    p = plot(
        dx,
        errors;
        label = "Combined ($(c.model))",
        color = CONDITION_COLORS[3],
        marker = :circle,
        linewidth = 2,
        xscale = :log10,
        yscale = :log10,
        xlabel = "dx (dimensionless)",
        ylabel = "Maximum absolute error",
        title = "Smooth analytic mode, t = 1",
        size = (800, 500),
        legend = :topleft,
    )
    plot!(
        p,
        dx,
        errors[1] .* dx ./ dx[1];
        label = "First order",
        color = :black,
        linestyle = :dot,
    )
    savefig(p, joinpath(directory, "convergence.png"))
end
const OUTPUT_NAMES =
    ("comparison.png", "conservation.png", "convergence.png", "summary.toml")

"""
    check_output_sizes(staged, output_dir, model)

置換される公式出力と保持する成果物を合わせて容量を検査する。

# 引数

- `staged`: 生成済み成果物の一時ディレクトリ。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `model`: `:linear` または `:nonlinear` の選択モデル。

# 返り値

成功時は `true`。1ファイル5MiB、N04全体10MiB、全課題100MiBを超えると停止する。
"""
function check_output_sizes(staged, output_dir, model)
    target = abspath(joinpath(output_dir, string(model)))
    replaced = Set(joinpath(target, n) for n in OUTPUT_NAMES)
    sizes = [filesize(joinpath(staged, n)) for n in OUTPUT_NAMES]
    all(s->0 < s <= 5 * 1024 ^ 2, sizes) || error("N04の1ファイル上限5 MiBを超えています")
    retained = 0
    if isdir(output_dir)
        for (root, _, files) in walkdir(output_dir), file in files
            path = abspath(joinpath(root, file))
            path in replaced && continue
            size = filesize(path)
            size <= 5 * 1024 ^ 2 || error("1ファイル上限5 MiBを超えています: $path")
            retained+=size
        end
    end
    sum(sizes) + retained <= 10 * 1024 ^ 2 || error("N04全体の上限10 MiBを超えています")
    # Standard course location: account for other assignments' existing results.
    exercises = dirname(@__DIR__)
    total = sum(sizes) + retained
    if basename(exercises) == "exercises" && isdir(exercises)
        for task in readdir(exercises; join = true)
            results = joinpath(task, "results")
            task == (@__DIR__) && continue
            isdir(results) || continue
            for (root, _, files) in walkdir(results), file in files
                total+=filesize(joinpath(root, file))
            end
        end
    end
    total <= 100 * 1024 ^ 2 || error("全課題の出力上限100 MiBを超えています")
end

"""
    install_output(source, destination)

生成済みファイルを指定された公式パスへコピーする。

# 引数

- `source`: コピー元のパス。
- `destination`: コピー先のパス。

# 返り値

コピー先のパス。同名ファイルを上書きする。
"""
install_output(source, destination) = cp(source, destination; force = true)

"""
    save_staged(draw, output_dir, model, summary)

一時領域で生成・検査した公式4ファイルを復元可能な形で反映する。

# 引数

- `draw`: 一時ディレクトリへ図を書き出す関数。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `model`: `:linear` または `:nonlinear` の選択モデル。
- `summary`: TOMLへ保存する診断情報。

# 返り値

`nothing`。生成・検査失敗は既存出力を保持し、反映失敗は触れたファイルを復元する。復元にも失敗するとバックアップの場所と両エラーを示し、バックアップを残す。
"""
function save_staged(draw, output_dir, model, summary)
    temporary = mktempdir(; cleanup = false)
    preserve_backup = false
    try
        staged = joinpath(temporary, "new")
        backup = joinpath(temporary, "backup")
        mkpath(staged)
        mkpath(backup)
        open(joinpath(staged, "summary.toml"), "w") do io
            TOML.print(io, summary; sorted = true)
        end
        draw(staged)
        check_output_sizes(staged, output_dir, model)
        for n in OUTPUT_NAMES[1:3]
            open(joinpath(staged, n)) do io
                read(io, 8) == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a] ||
                    error("PNG出力が不正です: $n")
            end
        end
        target = joinpath(output_dir, string(model))
        existing = Set{String}()
        for n in OUTPUT_NAMES
            path = joinpath(target, n)
            ispath(path) &&
                (!isfile(path) || islink(path)) &&
                error("公式出力先は通常ファイルである必要があります: $path")
            if isfile(path)
                cp(path, joinpath(backup, n))
                push!(existing, n)
            end
        end
        mkpath(target)
        touched = String[]
        try
            for n in OUTPUT_NAMES
                push!(touched, n)
                install_output(joinpath(staged, n), joinpath(target, n))
            end
        catch original
            try
                for n in touched
                    destination = joinpath(target, n)
                    if n in existing
                        cp(joinpath(backup, n), destination; force = true)
                    elseif isfile(destination)
                        rm(destination)
                    end
                end
            catch recovery
                preserve_backup = true
                error(
                    "出力反映と復元に失敗しました。バックアップ: $backup; 反映: $(sprint(showerror,original)); 復元: $(sprint(showerror,recovery))",
                )
            end
            rethrow(original)
        end
    finally
        preserve_backup || rm(temporary; recursive = true)
    end
end

"""
    write_outputs(output_dir, a, b, c, convergence)

比較結果と収束結果を診断TOMLと図へまとめる。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `a`: 移流のみの計算結果。
- `b`: 拡散のみの計算結果。
- `c`: 移流と拡散を含む計算結果。
- `convergence`: 滑らかな解析モードの格子収束結果。

# 返り値

`nothing`。選択モデルの公式出力だけを復元付きで反映し、無関係な追加ファイルを保持する。
"""
function write_outputs(output_dir, a, b, c, convergence)
    units =
        c.model == :linear ?
        "dimensionless temperature; integral is transported temperature content" :
        "dimensionless Burgers velocity; integral is velocity integral, not heat"

    # 計算値と条件を、保存する診断情報にまとめる。
    summary = Dict(
        "course_id" => "N04",
        "model" => string(c.model),
        "boundary" => "periodic",
        "domain" => [0., 2.],
        "units" => units,
        "initial" => "pulse",
        "t_final" => c.t_final,
        "requested_safety" => c.requested_safety,
        "advection_only" => summary_section(a),
        "diffusion_only" => summary_section(b),
        "combined" => summary_section(c),
        "convergence" => convergence,
    )
    save_staged(output_dir, c.model, summary) do directory
        make_plots(directory, a, b, c, convergence)
    end
end
