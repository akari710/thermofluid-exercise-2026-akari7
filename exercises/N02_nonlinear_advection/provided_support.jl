# 入力検証・描画・保存を担当する提供ファイル。学生の編集対象ではありません。
using Plots
using TOML
const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "results")

"""
    validate_boundary(boundary)

固定境界または周期境界の識別子を検証する。

# 引数

- `boundary`: 境界条件の識別子。

`boundary` は `:fixed` または `:periodic`。不正入力は `ArgumentError`。

# 返り値

検証に成功すると `true`。入力は変更しない。
"""
function validate_boundary(boundary)
    boundary in (:fixed, :periodic) ||
        throw(ArgumentError("boundaryは:fixedまたは:periodicです"))
end

"""
    validate_simulation_inputs(boundary, nx, cfl, t_final)

格子点数・CFL・最終時刻を検証する。

# 引数

- `boundary`: 境界条件の識別子。
- `nx`: x方向の格子点数。
- `cfl`: 指定するCFL数。
- `t_final`: 計算の最終時刻。

`nx` はBoolを除く3以上の整数、CFLと最終時刻は有限正値。不正入力は `ArgumentError`。

# 返り値

検証に成功すると `true`。入力は変更しない。
"""
function validate_simulation_inputs(boundary, nx, cfl, t_final)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_boundary(boundary)
    !(nx isa Bool) && nx >= 3 || throw(ArgumentError("nxは3以上の整数です（Bool不可）"))
    isfinite(cfl) && cfl > 0 || throw(ArgumentError("cflは有限な正の値です"))
    isfinite(t_final) && t_final > 0 || throw(ArgumentError("t_finalは有限な正の値です"))
end

"""
    validate_step_inputs(u_new, u_old, dt, dx, boundary)

浮動小数バッファの独立性・有限値・正の刻みを検証する。

# 引数

- `u_new`: 更新結果を書き込む配列。旧配列と重ならない独立したバッファ。
- `u_old`: 更新前の値を読む配列。更新中は変更しない。
- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。
- `dx`: x方向の有限な正の格子幅。
- `boundary`: 境界条件の識別子。

新旧は1始まり・同長・3点以上の浮動小数ベクトルで、互いに重ならない。
旧場は有限値、刻みと格子幅は有限正値。不正入力は書込み前に `ArgumentError`。

# 返り値

検証に成功すると `true`。旧配列を変更しない。CFL超過や負値自体は拒否しない。
"""
function validate_step_inputs(u_new, u_old, dt, dx, boundary)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_boundary(boundary)
    all(u -> u isa AbstractVector{<:AbstractFloat}, (u_new, u_old)) ||
        throw(ArgumentError("新旧バッファは浮動小数の一次元配列です"))
    Base.require_one_based_indexing(u_new, u_old)
    length(u_new) == length(u_old) >= 3 ||
        throw(ArgumentError("新旧は同長で3点以上必要です"))
    Base.mightalias(u_new, u_old) && throw(ArgumentError("新旧バッファが重なっています"))
    all(isfinite, u_old) || throw(ArgumentError("旧配列は有限値にしてください"))
    all(v -> isfinite(v) && v > 0, (dt, dx)) ||
        throw(ArgumentError("dtとdxは有限な正の値です"))
    # CFL超過や不安定化による負値でも，同じ後退差分を続けて破綻を観察する。
end

"""
    summary_section(result; periodic = false)

1境界条件の有限な診断量を辞書にまとめる。

# 引数

- `result`: シミュレーションの結果と診断量。
- `periodic`: 周期境界の計算結果、または周期表示を選ぶフラグ。

非有限な診断量は保存せずエラーにする。

# 返り値

刻み・CFL・極値・超過量の辞書。周期の場合は初期/最終の和とその変化も含む。
"""
function summary_section(result; periodic = false)
    # 計算値と条件を、保存する診断情報にまとめる。
    section = Dict{String,Any}("nx" => length(result.x))
    for key in (
        :dx,
        :steps,
        :dt,
        :max_cfl,
        :initial_minimum,
        :initial_maximum,
        :minimum,
        :maximum,
        :overshoot,
        :undershoot,
    )
        section[string(key)] = getproperty(result, key)
    end
    if periodic
        for key in (:initial_sum, :final_sum, :sum_change)
            section[string(key)] = getproperty(result, key)
        end
    end
    all(isfinite, values(section)) || error("非有限な診断量は保存できません")
    section
end

"""
    solution_plot(r; periodic = false)

初期値と最終値を比較する図を作る。

# 引数

- `r`: シミュレーションの結果と診断量。
- `periodic`: 周期境界の計算結果、または周期表示を選ぶフラグ。

# 返り値

Plotsの図。周期終点は表示専用に追加し、診断配列は変更しない。
"""
function solution_plot(r; periodic = false)
    label = periodic ? "Periodic boundary" : "Fixed boundary"
    # 周期の末尾へ表示用の点だけ追加する。診断の配列には混ぜない。
    x = periodic ? vcat(r.x, 2.0) : r.x
    initial = periodic ? vcat(r.u0, r.u0[1]) : r.u0
    final = periodic ? vcat(r.u, r.u[1]) : r.u
    low = min(minimum(r.u0), minimum(r.u))
    high = max(maximum(r.u0), maximum(r.u))
    padding = 0.05 * (high - low)
    limits = (min(0.9, low - padding), max(2.1, high + padding))
    p = plot(
        x,
        initial;
        label = "Initial",
        linewidth = 2,
        xlabel = "x",
        ylabel = "u",
        title = "$label, t = $(r.t_final)",
        size = (800, 500),
        ylims = limits,
    )
    plot!(p, x, final; label = "Final", linewidth = 2)
    return p
end

"""
    make_plots(output_dir, fixed, periodic)

固定境界と周期境界の公式図を保存する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `fixed`: 固定境界の計算結果。
- `periodic`: 周期境界の計算結果、または周期表示を選ぶフラグ。

# 返り値

`nothing`。指定ディレクトリへ `fixed-boundary.png` と `periodic.png` を保存する。
"""
function make_plots(output_dir, fixed, periodic)
    for (r, name, isperiodic) in
        ((fixed, "fixed-boundary.png", false), (periodic, "periodic.png", true))
        savefig(solution_plot(r; periodic = isperiodic), joinpath(output_dir, name))
    end
end

"""
    write_outputs(output_dir, fixed, periodic, cfl, t_final)

図と診断を一時領域で生成・検査してから公式パスへ反映する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `fixed`: 固定境界の計算結果。
- `periodic`: 周期境界の計算結果、または周期表示を選ぶフラグ。
- `cfl`: 指定するCFL数。
- `t_final`: 計算の最終時刻。

# 返り値

`nothing`。図2枚と `summary.toml` を保存する。生成・サイズ検査の失敗時は既存出力を保持する。反映中のコピー失敗には復元処理がない。
"""
function write_outputs(output_dir, fixed, periodic, cfl, t_final)
    summary = Dict(
        "course_id" => "N02",
        "domain" => [0.0, 2.0],
        "requested_cfl" => cfl,
        "t_final" => t_final,
        "fixed" => summary_section(fixed),
        "periodic" => summary_section(periodic; periodic = true),
    )
    # 全ファイルを一時領域で生成してから正式パスへ移す。生成失敗なら既存結果を保持。
    mktempdir() do temporary
        open(joinpath(temporary, "summary.toml"), "w") do io
            TOML.print(io, summary; sorted = true)
        end
        make_plots(temporary, fixed, periodic)
        names = ("fixed-boundary.png", "periodic.png", "summary.toml")
        sizes = [filesize(joinpath(temporary, n)) for n in names]
        maximum(sizes) <= 5 * 1024 ^ 2 && sum(sizes) <= 10 * 1024 ^ 2 ||
            error("N02出力サイズの上限を超えています")
        mkpath(output_dir)
        for name in names
            cp(joinpath(temporary, name), joinpath(output_dir, name); force = true)
        end
    end
end
