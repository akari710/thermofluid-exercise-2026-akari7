# N08/N09: 節点問題の生成 → ソルバ → 収束確認 → HDF5保存。

module N08N09Simulation
using ThermofluidExercise
include("provided_support.jl")
const E = ThermofluidExercise.Elliptic

"""
    simulate_case(
        id,
        nx,
        ny;
        laplace_solver = E.solve_laplace,
        poisson_solver = E.solve_poisson,
        atol = 1e-10,
        rtol = 1e-12,
        maxiter = 200000,
    )

独立に生成した問題を解き、収束済みの場と反復履歴をまとめる。

# 引数

- `id`: 格子点数と対応するケース名。
- `nx`: x方向の格子点数。
- `ny`: y方向の格子点数。
- `laplace_solver`: Laplace問題の収束情報と場を返す関数。
- `poisson_solver`: Poisson問題の収束情報と場を返す関数。
- `atol`: 有限で非負の絶対許容値。
- `rtol`: 有限で非負の相対許容値。
- `maxiter`: Boolを除く正整数の反復上限。

# 返り値

場・右辺・境界・座標・収束情報と反復履歴の辞書。未収束はID・格子・停止理由を示して拒否する。
"""
function simulate_case(
    id,
    nx,
    ny;
    laplace_solver = E.solve_laplace,
    poisson_solver = E.solve_poisson,
    atol = 1e-10,
    rtol = 1e-12,
    maxiter = 200000,
)
    p = problem(id, nx, ny)
    try
        result =
            id == "N08" ? laplace_solver(p.u0, p.g, p.dx, p.dy; atol, rtol, maxiter) :
            poisson_solver(p.u0, p.g, p.f, p.dx, p.dy; atol, rtol, maxiter)
        require(
            result.converged,
            "reason=$(result.reason), iterations=$(result.iterations), residual=$(last(result.residual_history)), threshold=$(result.threshold)",
        )
        Dict{String,Any}(
            "nx" => nx,
            "ny" => ny,
            "dx" => p.dx,
            "dy" => p.dy,
            "case_id" => p.case_id,
            "x" => p.x,
            "y" => p.y,
            "u" => result.u,
            "f" => p.f,
            "g" => p.g,
            "converged" => result.converged,
            "reason" => string(result.reason),
            "iterations" => result.iterations,
            "threshold" => result.threshold,
            "initial_residual" => first(result.residual_history),
            "residual_history" => result.residual_history,
            "update_history" => result.update_history,
            "iteration" => collect(0:result.iterations),
        )
    catch e
        error("$id $(p.case_id) $(nx)x$(ny): $(sprint(showerror,e))")
    end
end

"""
    main(;
        selection = "all",
        output_dir = DEFAULT_OUTPUT_DIR,
        laplace_solver = E.solve_laplace,
        poisson_solver = E.solve_poisson,
        publish_options...,
    )

選択IDの公式3格子を計算し、独立診断後に保存する。

# 引数

- `selection`: 処理する内容ID（N08、N09、all）。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `laplace_solver`: Laplace問題の収束情報と場を返す関数。
- `poisson_solver`: Poisson問題の収束情報と場を返す関数。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

nothing。選択IDのfields.h5だけを同じrun_idで復元付き反映する。未収束を保存しない。
"""
function main(;
    selection = "all",
    output_dir = DEFAULT_OUTPUT_DIR,
    laplace_solver = E.solve_laplace,
    poisson_solver = E.solve_poisson,
    publish_options...,
)
    ids = selection_ids(selection)
    run_id = bytes2hex(rand(UInt8, 16))
    names = Tuple(joinpath(id, "fields.h5") for id in ids)

    # 一時領域で処理を完了してから、公式出力を反映する。
    staged(output_dir, names; publish_options...) do stage
        for id in ids
            mkpath(joinpath(stage, id))
            cases = Dict(
                case_id(n...) =>
                    simulate_case(id, n...; laplace_solver, poisson_solver) for
                n in OFFICIAL_GRIDS
            )
            write_fields(
                joinpath(stage, id, "fields.h5"),
                cases,
                source_metadata(id, run_id),
            )
            diagnostics(
                read_fields(joinpath(stage, id, "fields.h5")),
                file_sha(joinpath(stage, id, "fields.h5")),
            )
        end
    end
    println(
        "$(join(ids,"・")) 定常場を保存しました。analyze.jl、plot.jlを実行してください。",
    )
end
abspath(PROGRAM_FILE) == (@__FILE__) && main(; selection = cli_selection(ARGS))
end
