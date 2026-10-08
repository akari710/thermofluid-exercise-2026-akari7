# N08/N09: 保存した発展場・残差履歴・格子収束を描画します。

using Plots

"""
    plot_extension(dir, pair, summary)

発展の保存場・残差履歴・格子収束を作図する。

# 引数

- `dir`: 保存場を読むディレクトリ。
- `pair`: 同じrun_idとコード来歴を持つLaplace/Poissonの保存場。発展課題では選択した発展の保存場。
- `summary`: TOMLへ保存する診断情報。

# 返り値

表示設定とrun_id・保存場/解析/PNGのhashを持つ辞書。dirに公式3図を保存する。
"""
function plot_extension(dir, pair, summary)
    ext = pair.metadata["extension_id"]
    ids = ext == "neumann" ? collect(NEUMANN_CASES) : ["jacobi", "gauss_seidel", "sor"]
    panels = Any[]
    config = Dict{String,Any}()
    for id in ids
        c = pair.cases[id * "_" * case_id(33, 25)]
        p = ext == "neumann" ? neumann_problem(id, 33, 25) : problem("N09", 33, 25)
        exact = copy(p.exact)
        ext == "neumann" && is_pure(p.bc) && (exact.-=weighted_mean(exact, p.dx, p.dy))
        limits = extrema(vcat(vec(c["u"]), vec(exact)))
        error = c["u"] - exact
        e = max(maximum(abs, error), eps())
        # GR uses its own automatic colorbar ticks. Display a stated power-of-ten
        # unit so tiny errors retain readable ticks and a symmetric physical range.
        exponent = floor(Int, log10(e))
        scale = 10.0 ^ exponent
        config[id] = Dict(
            "field_color_range" => collect(limits),
            "error_color_range" => [-e, e],
            "error_display_scale" => scale,
            "error_display_color_range" => [-e / scale, e / scale],
        )
        for (a, t, clims, color) in (
            (c["u"], "$id numerical", limits, :viridis),
            (exact, "Exact", limits, :viridis),
            (error ./ scale, "Error (×10^$exponent)", (-e / scale, e / scale), :RdBu),
        )
            push!(
                panels,
                heatmap(
                    p.x,
                    p.y,
                    transpose(a);
                    title = t,
                    xlabel = "x",
                    ylabel = "y",
                    aspect_ratio = :equal,
                    clims,
                    color,
                    titlefontsize = 9,
                    tickfontsize = 6,
                ),
            )
        end
    end
    savefig(
        plot(
            panels...;
            layout = (length(ids), 3),
            size = (1200, 350length(ids)),
            margin = 5Plots.mm,
        ),
        joinpath(dir, "fields.png"),
    )
    residual = plot(;
        yscale = :log10,
        xlabel = "Iteration",
        ylabel = "Linf residual / update",
        size = (1000, 600),
        legend = :outerright,
    )
    for (id, color) in zip(ids, (:blue, :orange, :green, :purple, :red))
        c = pair.cases[id * "_" * case_id(33, 25)]
        floor = 1e-16
        plot!(
            residual,
            c["iteration"],
            max.(c["residual_history"], floor);
            label = id,
            color,
            linewidth = 2,
        )
        plot!(
            residual,
            c["iteration"],
            max.(c["update_history"], floor);
            label = "$id update",
            color,
            linestyle = :dash,
        )
        plot!(
            residual,
            [0, c["iterations"]],
            fill(max(c["threshold"], floor), 2);
            label = "$id threshold",
            color,
            linestyle = :dot,
        )
    end
    savefig(residual, joinpath(dir, "residual.png"))
    convergence = plot(;
        xscale = :log10,
        yscale = :log10,
        xlabel = "dx",
        ylabel = "RMS error",
        size = (750, 500),
    )
    order_ids = ext == "neumann" ? ["pure_poisson_zero_flux"] : ids
    hs = [2 / (n[1] - 1) for n in OFFICIAL_GRIDS]
    for id in order_ids
        errors =
            [summary["cases"][id * "_" * case_id(n...)]["l2_error"] for n in OFFICIAL_GRIDS]
        plot!(convergence, hs, errors; marker = :circle, label = id)
        plot!(
            convergence,
            hs,
            errors[1] .* (hs ./ hs[1]) .^ 2;
            label = "$id second order",
            linestyle = :dash,
        )
    end
    savefig(convergence, joinpath(dir, "convergence.png"))
    merge(
        config,
        Dict(
            "display_grid" => case_id(33, 25),
            "axis_order" => "y,x",
            "log_floor" => 1e-16,
            "schema_version" => 1,
            "extension_id" => ext,
            "task_id" => "N09",
            "run_id" => summary["run_id"],
            "source_fields_sha256" => file_sha(joinpath(dir, "fields.h5")),
            "source_summary_sha256" => file_sha(joinpath(dir, "summary.toml")),
            "png_sha256" => Dict(n => file_sha(joinpath(dir, n)) for n in PNG_FILES),
        ),
    )
end
