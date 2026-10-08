# 入力検証・描画・保存の提供ファイル。学生の編集対象ではありません。
# run.jl から読み込まれ、N03Diffusion モジュールの一部になります。
# 関数は「入力検証 → 計算結果の整理 → 描画 → 保存」の順に並んでいます。

using Plots, TOML

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "results")

"""
    positive_finite(v, name)

指定された値が有限な正の値であることを確認する。

# 引数

- `v`: 検査する実数値。
- `name`: エラーメッセージに表示する値の名前。

非有限値または0以下の値なら、`name` を含む `ArgumentError` を投げる。

# 返り値

条件を満たす場合は `nothing`。
"""
positive_finite(v, name) = isfinite(v) && v > 0 ? nothing : throw(ArgumentError("$name は有限な正の値です"))

"""
    validate_boundary(b)

境界条件が `:fixed` または `:insulated` であることを確認する。

# 引数

- `b`: 検査する境界条件の指定。それ以外の指定では `ArgumentError` を投げる。

# 返り値

条件を満たす場合は `nothing`。
"""
validate_boundary(b) = b in (:fixed, :insulated) ? nothing : throw(ArgumentError("boundaryは:fixedまたは:insulatedです"))

"""
    validate_initial(i)

初期条件が `:pulse` または `:mode` であることを確認する。

# 引数

- `i`: 検査する初期条件の指定。それ以外の指定では `ArgumentError` を投げる。

# 返り値

条件を満たす場合は `nothing`。
"""
validate_initial(i) = i in (:pulse, :mode) ? nothing : throw(ArgumentError("initialは:pulseまたは:modeです"))

"""
    validate_buffers(a, b, boundary)

新旧の温度バッファを安全に更新できることを確認する。

# 引数

- `a`: 更新先の浮動小数ベクトル。この検証関数では変更しない。
- `b`: 旧状態の浮動小数ベクトル。全要素を有限値とし、この検証関数では変更しない。
- `boundary`: 固定温度の `:fixed`、または断熱の `:insulated`。

新旧配列は1始まり・同長・3点以上とし、互いに重ならないようにする。
条件を満たさなければ `ArgumentError` を投げる。

# 返り値

条件を満たす場合は `true`。呼出し元は返り値を使わず、検査を通過したことを利用する。
"""
function validate_buffers(a, b, boundary)
    # 境界条件と配列の形を、値を書き換える前に検査する。
    validate_boundary(boundary)
    all(u -> u isa AbstractVector{<:AbstractFloat}, (a, b)) ||
        throw(ArgumentError("新旧バッファは浮動小数の一次元配列です"))
    Base.require_one_based_indexing(a, b)
    length(a) == length(b) >= 3 || throw(ArgumentError("新旧は同長で3点以上必要です"))

    # 新旧の領域が重なると、更新済みの値を旧状態として読んでしまう。
    Base.mightalias(a, b) && throw(ArgumentError("新旧バッファが重なっています"))
    all(isfinite, b) || throw(ArgumentError("旧配列は有限値にしてください"))
end

"""
    validate_simulation_inputs(boundary, nx, d, fo, t, initial)

シミュレーションの境界・初期条件、格子点数、時間刻みに関係する値を検査する。

# 引数

- `boundary`: 固定温度の `:fixed`、または断熱の `:insulated`。
- `nx`: 格子点数。呼出し元で整数に制限され、ここでは3以上かつ `Bool` 以外を求める。
- `d`: 有限な正の拡散係数。
- `fo`: 有限な正の要求Fourier数。
- `t`: 有限な正の終了時刻。
- `initial`: 矩形分布の `:pulse`、または解析モードの `:mode`。

条件を満たさなければ `ArgumentError` を投げる。

# 返り値

条件を満たす場合は `nothing`。
"""
function validate_simulation_inputs(boundary, nx, d, fo, t, initial)
    # 境界条件・初期条件と格子点数を検査する。
    validate_boundary(boundary)
    validate_initial(initial)
    !(nx isa Bool) && nx >= 3 || throw(ArgumentError("nxは3以上の整数です（Bool不可）"))

    # 時間刻みの決定に必要な値を検査する。
    for (v, name) in ((d, "diffusivity"), (fo, "fo"), (t, "t_final"))
        positive_finite(v, name)
    end
end

"""
    ensure_finite(values, step, t, fo)

温度や診断量の非有限値を検出し、計算状況を示して停止する。

# 引数

- `values`: 有限性を確認する値の配列またはタプル。
- `step`: エラーメッセージに表示するステップ番号。
- `t`: エラーメッセージに表示する時刻。
- `fo`: エラーメッセージに表示する実効Fourier数。

`values` に非有限値があれば、ステップ番号・時刻・実効Fourier数を含むエラーで停止する。

# 返り値

全ての値が有限なら `true`。呼出し元は返り値を使わず、検査を通過したことを利用する。
"""
function ensure_finite(values, step, t, fo)
    all(isfinite, values) || error("非有限値のため停止: step=$step, t=$t, effective_fo=$fo")
end

"""
    summary_section(r)

シミュレーション結果から、TOML保存用の診断量を取り出す。

# 引数

- `r`: `simulate` が返す計算結果の `NamedTuple`。

保存対象の診断量に非有限値があれば、保存前にエラーで停止する。

# 返り値

格子点数・計算条件・熱量・温度の極値を文字列キーでまとめた `Dict{String,Any}`。
温度配列や時刻履歴は含めない。
"""
function summary_section(r)
    # 配列全体ではなく、計算条件と診断量を保存する。
    section = Dict{String,Any}("nx" => length(r.x))
    for key in (
        :dx, :dt, :steps, :requested_fo, :fo, :t_final, :diffusivity,
        :initial_heat, :final_heat, :heat_change,
        :initial_minimum, :initial_maximum, :minimum, :maximum,
    )
        section[string(key)] = getproperty(r, key)
    end

    all(isfinite, values(section)) || error("非有限な診断量は保存できません")
    return section
end

"""
    convergence_results()

各境界条件の解析解モードを3段階の格子で計算し、格子収束を調べる。

# 引数

なし。41・81・161点の格子を使い、その他の計算条件は `simulate` の既定値に従う。
実行にはN03のTODO 3か所の実装が必要。
誤差または観測次数が非有限値になった場合はエラーで停止する。

# 返り値

`"fixed"` と `"insulated"` をキーに持つ `Dict{String,Any}`。
各境界条件の辞書には、格子点数 `"nx"`、格子幅 `"dx"`、最大絶対誤差 `"errors"`、観測次数 `"orders"` を含む。
"""
function convergence_results()
    result = Dict{String,Any}()
    for boundary in (:fixed, :insulated)
        # 格子を細かくして、同じ終了時刻の数値解と解析解を比較する。
        nx = [41, 81, 161]
        runs = [simulate(; boundary, nx=n, initial=:mode) for n in nx]
        errors = [maximum(abs.(r.u - analytic_solution(r.x, r.t_final; boundary))) for r in runs]

        # 隣り合う格子の誤差比から観測次数を求める。
        orders = log2.(errors[1:2] ./ errors[2:3])
        all(isfinite, vcat(errors, orders)) || error("非有限な収束診断量です")
        result[string(boundary)] = Dict(
            "nx" => nx,
            "dx" => [r.dx for r in runs],
            "errors" => errors,
            "orders" => orders,
        )
    end

    return result
end

# 比較図では時刻を色で、境界条件を線種で区別する。
const SNAPSHOT_TIMES = (0.25, 0.50, 0.75, 1.0)
# 白い背景でも見えるよう、viridis の淡い黄色の端を避ける。
const SNAPSHOT_COLORS = [get(cgrad(:viridis), v) for v in (0.08, 0.30, 0.53, 0.76)]

"""
    boundary_snapshots(final; boundary)

標準パルスの比較図に使う、4つの時刻の計算結果を用意する。

# 引数

- `final`: 標準計算の最終結果。同じ時刻の計算結果には、この値を再利用する。
- `boundary`: 固定温度の `:fixed`、または断熱の `:insulated`。必須で指定する。

格子・初期値・時間刻み・実効Fourier数が `final` と一致しなければエラーで停止する。
温度に非有限値がある場合も停止する。

# 返り値

時刻0.25・0.50・0.75・1.0の計算結果を順に持つ4要素のタプル。
各要素は `simulate` と同じ形式の `NamedTuple`。
"""
function boundary_snapshots(final; boundary)
    return map(SNAPSHOT_TIMES) do t
        # 比較する各時刻で、最終計算と同じ条件を使う。
        r = t == final.t_final ? final : simulate(
            ; boundary, nx=length(final.x), diffusivity=final.diffusivity,
            fo=final.requested_fo, t_final=t,
        )
        r.x == final.x && r.u0 == final.u0 && r.dt == final.dt && r.fo == final.fo ||
            error("比較図の時刻・格子・刻みが標準計算と一致しません: t=$t")
        ensure_finite(r.u, r.steps, r.t_final, r.fo)
        r
    end
end

"""
    boundary_comparison_plot(fixed, insulated)

固定温度と断熱の温度分布を、時刻の色と境界条件の線種で比較する図を返す。

# 引数

- `fixed`: 固定温度の標準計算で `simulate` が返した結果。
- `insulated`: 断熱の標準計算で `simulate` が返した結果。

# 返り値

温度分布のパネルと、時刻・線種の凡例を組み合わせた `Plots.Plot`。
ファイル保存は呼出し元で行う。
"""
function boundary_comparison_plot(fixed, insulated)
    # 温度分布のパネルに、初期値と各時刻の数値解を重ねる。
    p = plot(
        fixed.x, fixed.u0;
        label="", color=:gray, linestyle=:dot, linewidth=2,
        xlabel="x (dimensionless)", ylabel="u (dimensionless)",
        title="Diffusion: fixed / insulated", ylims=(-0.05, 1.05), legend=false,
    )
    for (final, boundary, style) in ((fixed, :fixed, :solid), (insulated, :insulated, :dash))
        for (r, color) in zip(boundary_snapshots(final; boundary), SNAPSHOT_COLORS)
            plot!(p, r.x, r.u; label="", color, linestyle=style, linewidth=2)
        end
    end

    # 凡例を図の右側に分けて置き、曲線を隠さずに時刻を示す。
    time_key = plot(; axis=false, grid=false, legend=:left, legendtitle="Time", framestyle=:none)
    for (label, color) in zip(("t = 0.25", "t = 0.50", "t = 0.75", "t = 1.0"), SNAPSHOT_COLORS)
        scatter!(time_key, [NaN], [NaN]; label, color, markershape=:square, markerstrokewidth=0)
    end

    # 線種の見本を長く描き、破線と点線を見分けやすくする。
    boundary_key = plot(
        ; axis=false, grid=false, legend=false, framestyle=:none,
        xlims=(0, 1), ylims=(0, 1),
    )
    annotate!(boundary_key, 0.02, 0.85, text("Line style", 12, :left))
    for (y, label, style, color) in (
        (0.68, "Fixed temperature", :solid, :black),
        (0.50, "Insulated", :dash, :black),
        (0.32, "Initial (t = 0)", :dot, :gray),
    )
        plot!(boundary_key, [0.02, 0.28], [y, y]; label="", color, linestyle=style, linewidth=2)
        annotate!(boundary_key, 0.33, y, text(label, 9, :left))
    end

    return plot(
        p, time_key, boundary_key;
        layout=@layout([a{0.72w} [b; c]]), size=(800, 500), background_color=:white,
    )
end

"""
    comparison_plot(a, b; unstable=false)

初期分布と2つの計算結果を重ねた図を返す。

# 引数

- `a`: 初期分布と1つ目の最終分布に使う `simulate` の結果。
- `b`: 2つ目の最終分布に使う `simulate` の結果。`a` と同じ終了時刻で比較する。
- `unstable`: `true` なら安定・不安定の区別と実効Fourier数を凡例に表示する。
  `false` なら固定温度・断熱の比較として表示する。

# 返り値

初期分布と2つの最終分布を重ねた `Plots.Plot`。
縦軸は全ての分布が収まる範囲に調整する。ファイル保存は呼出し元で行う。
"""
function comparison_plot(a, b; unstable=false)
    # 初期分布と最終分布を重ねる。
    p = plot(
        a.x, a.u0;
        label="Initial", linewidth=2, xlabel="x (dimensionless)", ylabel="u (dimensionless)",
        title="Diffusion, t = $(a.t_final)", size=(800, 500),
    )
    plot!(p, a.x, a.u; label=unstable ? "Stable (fo=$(round(a.fo; digits=4)))" : "Fixed temperature", linewidth=2)
    plot!(p, b.x, b.u; label=unstable ? "Unstable (fo=$(round(b.fo; digits=4)))" : "Insulated", linewidth=2)

    # 不安定計算で値が増幅しても、全ての分布が縦軸に収まるようにする。
    low = min(minimum(a.u0), minimum(a.u), minimum(b.u))
    high = max(maximum(a.u0), maximum(a.u), maximum(b.u))
    padding = max(0.05 * (high - low), 0.01)
    ylims!(p, (low - padding, high + padding))
    return p
end

"""
    make_plots(directory, fixed, insulated, convergence)

温度分布・熱量履歴・格子収束の3図を、`directory` にPNGで保存する。

# 引数

- `directory`: 図を保存する既存のディレクトリ。
- `fixed`: 固定温度の標準計算で `simulate` が返した結果。
- `insulated`: 断熱の標準計算で `simulate` が返した結果。
- `convergence`: `convergence_results` が返した格子収束の診断辞書。

# 返り値

最後に保存した `convergence.png` の絶対パス。
同じディレクトリに `boundary-comparison.png` と `heat-content.png` も保存する。
"""
function make_plots(directory, fixed, insulated, convergence)
    # 境界条件による温度分布の違いを保存する。
    savefig(boundary_comparison_plot(fixed, insulated), joinpath(directory, "boundary-comparison.png"))

    # 台形則で求めた熱量の時間変化を比較する。
    p = plot(
        fixed.times, fixed.heat_history;
        label="Fixed temperature", color=:black, linestyle=:solid, linewidth=2,
        xlabel="t (dimensionless)", ylabel="Heat content H (dimensionless)", size=(800, 500),
    )
    plot!(p, insulated.times, insulated.heat_history; label="Insulated", color=:black, linestyle=:dash, linewidth=2)
    hline!(p, [fixed.initial_heat]; label="Initial H", color=:gray, linestyle=:dot)
    savefig(p, joinpath(directory, "heat-content.png"))

    # 格子幅と誤差の関係を、両対数軸で比較する。
    p = plot(
        ; xscale=:log10, yscale=:log10, xlabel="dx", ylabel="Maximum absolute error",
        size=(800, 500), legend=:topleft,
    )
    for (boundary, label, style, marker) in (
        ("fixed", "Fixed temperature", :solid, :circle),
        ("insulated", "Insulated", :dash, :diamond),
    )
        c = convergence[boundary]
        plot!(p, c["dx"], c["errors"]; label, color=:black, linestyle=style, marker, linewidth=2)
    end

    # 固定温度の粗い格子の誤差を基準に、2次収束の参照線を添える。
    c = convergence["fixed"]
    plot!(p, c["dx"], c["errors"][1] .* (c["dx"] ./ c["dx"][1]).^2; label="Second order", color=:black, linestyle=:dot)
    savefig(p, joinpath(directory, "convergence.png"))
end

"""
    save_staged(draw, output_dir, summary, names)

TOMLと図を一時ディレクトリに生成し、サイズとPNG署名を検査して保存する。

# 引数

- `draw`: `draw(directory)` として呼び出す描画関数。渡された一時ディレクトリに図を生成する。
- `output_dir`: 検査後のファイルをコピーする出力先。
- `summary`: `summary.toml` に書き込む診断辞書。
- `names`: `summary.toml` を含む保存対象のファイル名の列。全て一時ディレクトリに生成する必要がある。

各ファイルは0より大きく5MiB以下、合計は10MiB以下とする。PNGファイルは署名も検査する。
生成・検査中の失敗では既存の出力を変更しない。
コピー中の失敗では、一部のファイルが更新済みになる場合がある。

# 返り値

保存が完了した場合は `nothing`。出力先がなければ作成し、指定されたファイルを順に上書きする。
"""
function save_staged(draw, output_dir, summary, names)
    mktempdir() do temporary
        # まず全ての出力を一時ディレクトリに生成する。
        open(joinpath(temporary, "summary.toml"), "w") do io
            TOML.print(io, summary; sorted=true)
        end
        draw(temporary)

        # 空ファイルやサイズ超過、不正なPNGを保存前に検出する。
        sizes = [filesize(joinpath(temporary, n)) for n in names]
        all(s -> 0 < s <= 5 * 1024^2, sizes) && sum(sizes) <= 10 * 1024^2 ||
            error("N03出力サイズの上限を超えています")
        for n in filter(n -> endswith(n, ".png"), names)
            open(joinpath(temporary, n)) do io
                read(io, 8) == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a] ||
                    error("PNG出力が不正です: $n")
            end
        end

        # 検査を通過したファイルを、指定された順に反映する。
        mkpath(output_dir)
        for n in names
            cp(joinpath(temporary, n), joinpath(output_dir, n); force=true)
        end
    end
end

"""
    write_outputs(output_dir, fixed, insulated, convergence)

標準計算の診断量を `summary.toml` にまとめ、3図とともに検査して保存する。

# 引数

- `output_dir`: 公式4出力の保存先。
- `fixed`: 固定温度の標準計算で `simulate` が返した結果。
- `insulated`: 断熱の標準計算で `simulate` が返した結果。
- `convergence`: `convergence_results` が返した格子収束の診断辞書。

保存時の失敗については [`save_staged`](@ref) を参照する。

# 返り値

保存が完了した場合は `nothing`。
出力先に3図のPNGと `summary.toml` を保存する。
"""
function write_outputs(output_dir, fixed, insulated, convergence)
    # TOMLの共通条件と、各境界条件の診断量をまとめる。
    summary = Dict(
        "course_id" => "N03",
        "units" => "dimensionless; heat_content H uses trapezoidal endpoint weights",
        "domain" => [0., 2.], "diffusivity" => 0.1, "requested_fo" => 0.4, "t_final" => 1.,
        "fixed" => summary_section(fixed),
        "insulated" => summary_section(insulated),
        "convergence" => convergence,
    )

    # doブロックで渡す描画関数が、save_staged の第1引数になる。
    save_staged(output_dir, summary, (
        "boundary-comparison.png", "heat-content.png", "convergence.png", "summary.toml",
    )) do temp
        make_plots(temp, fixed, insulated, convergence)
    end
end

"""
    make_stability_plot(directory, stable, unstable)

安定・不安定計算の比較図を、`directory` にPNGで保存する。

# 引数

- `directory`: 図を保存する既存のディレクトリ。
- `stable`: 安定条件で `simulate` が返した結果。
- `unstable`: 不安定条件で `simulate` が返した結果。`stable` と同じ終了時刻で比較する。

# 返り値

保存した `stability-comparison.png` の絶対パス。
"""
function make_stability_plot(directory, stable, unstable)
    savefig(comparison_plot(stable, unstable; unstable=true), joinpath(directory, "stability-comparison.png"))
end

"""
    write_stability_outputs(output_dir, stable, unstable)

安定・不安定計算の診断量と比較図を、検査して保存する。

# 引数

- `output_dir`: 任意実験の出力の保存先。
- `stable`: 安定条件で `simulate` が返した結果。
- `unstable`: 不安定条件で `simulate` が返した結果。`stable` と同じ終了時刻で比較する。

保存時の失敗については [`save_staged`](@ref) を参照する。

# 返り値

保存が完了した場合は `nothing`。
出力先に `stability-comparison.png` と `summary.toml` を保存する。
"""
function write_stability_outputs(output_dir, stable, unstable)
    # 2つの計算結果を、同じ固定温度条件の実験として記録する。
    summary = Dict(
        "course_id" => "N03", "units" => "dimensionless", "boundary" => "fixed",
        "stable" => summary_section(stable), "unstable" => summary_section(unstable),
    )

    save_staged(output_dir, summary, ("stability-comparison.png", "summary.toml")) do temp
        make_stability_plot(temp, stable, unstable)
    end
end
