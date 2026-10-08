# N08/N09: 混合境界 → 純Neumannの計算・独立診断・作図を発展出力へ保存。

module NeumannRun
include("neumann.jl")
include("extension_support.jl")
include("extension_plots.jl")

"""
    main(;
        output_dir = DEFAULT_OUTPUT_DIR,
        solver = NeumannExtension.solve_neumann,
        publish_options...,
    )

混合境界から純Neumannへ公式発展ケースを計算・検査・保存する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `solver`: 収束情報と場を返す発展問題のソルバ。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

nothing。N09/extensions/neumann/の公式6ファイルだけを復元付き反映する。適合条件違反・未収束・診断失敗では反映しない。
"""
function main(;
    output_dir = DEFAULT_OUTPUT_DIR,
    solver = NeumannExtension.solve_neumann,
    publish_options...,
)
    cases = Dict{String,Any}()
    run_id = bytes2hex(rand(UInt8, 16))
    # All compatibility checks and calculations precede staging and publication.
    for id in NEUMANN_CASES, (nx, ny) in OFFICIAL_GRIDS
        p = neumann_problem(id, nx, ny)
        result = solver(p.u0, p.f, p.dx, p.dy, p.bc)
        cases[id * "_" * case_id(nx, ny)] = extension_case(p, result, id)
    end
    relative = joinpath("N09", "extensions", "neumann")
    names = Tuple(joinpath(relative, n) for n in OUTPUT_FILES)

    # 一時領域で処理を完了してから、公式出力を反映する。
    staged(output_dir, names; publish_options...) do stage
        dir = joinpath(stage, relative)
        mkpath(dir)
        path = joinpath(dir, "fields.h5")
        pair = write_extension_fields(path, cases, extension_metadata("neumann", run_id))
        summary = extension_diagnostics(pair, file_sha(path))
        write_toml(joinpath(dir, "summary.toml"), summary)
        write_toml(joinpath(dir, "plots.toml"), plot_extension(dir, pair, summary))
        check_extension(dir)
    end
    println("Neumann発展の混合→純Neumannを検証・保存しました")
end
if abspath(PROGRAM_FILE) == (@__FILE__)
    isempty(ARGS) || error("run_neumann.jlは引数なしです")
    main()
end
end
