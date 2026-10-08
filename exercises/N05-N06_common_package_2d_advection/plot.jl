# N05/N06: 保存場と解析の出自確認 → 表示設定 → 図とplots.toml保存。表示設定は受講生が調整できます。

module N06Plots
include("provided_support.jl")
using Plots
# 学生の編集対象: 物理条件・診断値を変えずに表示だけを調整する。
const COLORMAP = :viridis
const COLOR_LIMITS = (0.4, 1.6)
const DISPLAY_CASE = "n080x060"

"""
    draw(
        stage,
        data,
        s;
        colormap = COLORMAP,
        color_limits = COLOR_LIMITS,
        display_case = DISPLAY_CASE,
    )

保存場と解析結果から分布・保存量・収束の図を描く。

# 引数

- `stage`: 新しい成果物を生成する一時ディレクトリ。
- `data`: 検証済みの保存場とmetadata。
- `s`: 保存場と出自を照合した解析結果。
- `colormap`: 場の表示に使うカラーマップ。
- `color_limits`: 場の表示範囲(lower, upper)。有限でlower<upper。
- `display_case`: 再描画する保存ケース名。

# 返り値

最後のsavefigの結果。一時出力先に図3枚を保存する。物理条件と保存場は変更しない。
"""
function draw(
    stage,
    data,
    s;
    colormap = COLORMAP,
    color_limits = COLOR_LIMITS,
    display_case = DISPLAY_CASE,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        length(color_limits) == 2 &&
            all(isfinite, color_limits) &&
            color_limits[1] < color_limits[2],
        "色範囲が不正です",
    )
    require(haskey(data.cases, display_case), "表示ケースがありません: $display_case")
    c = data.cases[display_case]
    time = last(c["time"])
    initial = c["u"][:, :, 1]
    final = c["u"][:, :, end]
    exact = exact_field(
        c["x"],
        c["y"],
        time;
        cx = data.metadata["cx"],
        cy = data.metadata["cy"],
    )
    error = final - exact
    panels = []
    for (u, title) in (
        (initial, "Initial, t=0"),
        (final, "Numerical, t=$time"),
        (exact, "Exact, t=$time"),
        (error, "Numerical - exact"),
    )
        iserror = u===error
        amplitude = max(maximum(abs, error), eps())
        push!(
            panels,
            heatmap(
                c["x"],
                c["y"],
                permutedims(u);
                xlabel = "x",
                ylabel = "y",
                title,
                color = iserror ? :RdBu : colormap,
                clims = iserror ? (-amplitude, amplitude) : color_limits,
                xlims = (0., 2.),
                ylims = (0., 1.),
                aspect_ratio = :equal,
            ),
        )
    end
    savefig(
        plot(panels...; layout = (2, 2), size = (1100, 660), margin = 5Plots.mm),
        joinpath(stage, "fields.png"),
    )
    mass =
        plot(; xlabel = "Saved time", ylabel = "Mass change / 1e-15", legend = :bottomleft)
    errors = plot(; xlabel = "Saved time", ylabel = "L2 error", legend = :topleft)
    for id in s["convergence"]["case_ids"]
        d = s["cases"][id]
        plot!(
            mass,
            d["time"],
            (d["mass"] .- first(d["mass"])) ./ 1e-15;
            label = id,
            marker = :circle,
        )
        plot!(errors, d["time"], d["l2_error"]; label = id, marker = :circle)
    end
    savefig(
        plot(mass, errors; layout = (1, 2), size = (1000, 400), margin = 5Plots.mm),
        joinpath(stage, "diagnostics.png"),
    )
    ids = s["convergence"]["case_ids"]
    h = [s["cases"][id]["dx"] for id in ids]
    err = s["convergence"]["errors"]
    p = plot(
        h,
        err;
        xscale = :log10,
        yscale = :log10,
        xlabel = "dx (dy refined together)",
        ylabel = "Final L2 error",
        label = "Numerical",
        marker = :circle,
        legend = :topleft,
        size = (700, 460),
    )
    plot!(p, h, first(err) .* h ./ first(h); label = "First order", linestyle = :dash)
    savefig(p, joinpath(stage, "convergence.png"))
end

"""
    main(;
        input_path = joinpath(DEFAULT_OUTPUT_DIR, FIELD_NAME),
        summary_path = joinpath(DEFAULT_OUTPUT_DIR, "summary.toml"),
        output_dir = DEFAULT_OUTPUT_DIR,
        colormap = COLORMAP,
        color_limits = COLOR_LIMITS,
        display_case = DISPLAY_CASE,
        publish_options...,
    )

保存場と解析のhashを照合し、完成した診断から再作図する。

# 引数

- `input_path`: 読取り専用の保存場HDF5のパス。
- `summary_path`: 保存場と対応する解析TOMLのパス。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `colormap`: 場の表示に使うカラーマップ。
- `color_limits`: 場の表示範囲(lower, upper)。有限でlower<upper。
- `display_case`: 再描画する保存ケース名。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

最後のprintlnの結果(nothing)。図3枚とplots.tomlを復元付きで反映する。分散未完成や古い解析結果は反映前に拒否する。
"""
function main(;
    input_path = joinpath(DEFAULT_OUTPUT_DIR, FIELD_NAME),
    summary_path = joinpath(DEFAULT_OUTPUT_DIR, "summary.toml"),
    output_dir = DEFAULT_OUTPUT_DIR,
    colormap = COLORMAP,
    color_limits = COLOR_LIMITS,
    display_case = DISPLAY_CASE,
    publish_options...,
)
    # 保存済みの値と出自を検証してから、解析や作図に使う。
    data = read_fields(input_path)
    s = read_summary(summary_path, input_path)

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(get(s, "diagnostics_complete", false), "spatial_varianceが未実装です")
    names = ("fields.png", "diagnostics.png", "convergence.png", "plots.toml")

    # 一時領域で処理を完了してから、公式出力を反映する。
    staged(output_dir, names; publish_options...) do stage
        draw(stage, data, s; colormap, color_limits, display_case)
        write_toml(
            joinpath(stage, "plots.toml"),
            Dict(
                "schema_version" => 1,
                "source_fields_sha256" => file_sha(input_path),
                "source_summary_sha256" => file_sha(summary_path),
                "display_case" => display_case,
                "time" => [
                    first(data.cases[display_case]["time"]),
                    last(data.cases[display_case]["time"]),
                ],
                "colormap" => string(colormap),
                "color_limits" => collect(color_limits),
                "error_colormap" => "RdBu",
                "mass_display_scale" => 1e-15,
            ),
        )
    end
    println("N06を保存場から再作図しました。")
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
