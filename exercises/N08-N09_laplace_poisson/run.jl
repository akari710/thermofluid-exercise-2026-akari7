# N08/N09の提供済み一括入口。必修の編集先はsrc/N08N09Elliptic.jlです。
# N08の4関数: apply_dirichlet! / laplace_jacobi_step! / laplace_residual! / residual_converged。
# N09の追加2関数: poisson_jacobi_step! / poisson_residual!。解析・保存・作図は提供済みです。
# tests.jlの配布必須テストを読み、同じファイルの「自作」枠に自分のテストを書きます。
# 読む順: 課題N08/N09 → 必修関数とtests.jl → 提供のsolve_driver → simulate.jl → analyze.jl → plot.jl → このmain。
# 任意の発展はextensions/neumann.jl / relaxation.jlと各発展テスト・実行入口を別に読みます。
# 計算 → 保存場の独立診断 → 作図 → 完了検査を一時領域で終え、選択IDの出力へ反映します。

module N08N09Run
include("simulate.jl")
include("analyze.jl")
include("plot.jl")

"""
    main(;
        selection = "all",
        output_dir = N08N09Simulation.DEFAULT_OUTPUT_DIR,
        simulation = N08N09Simulation.main,
        analysis = N08N09Analysis.main,
        plotting = N08N09Plots.main,
        publish_options...,
    )

一時領域で計算・診断・作図・完了検査を行い、公式出力を反映する。
既定の `selection="all"` には必修6関数、`selection="N08"` にはN08の4関数の実装が必要。
任意の `extensions/` はこの入口では実行しない。
TODOの未実装や反復の非収束では生成・保存が止まり、公式出力の反映へ進まない。

# 引数

- `selection`: 処理する内容ID（N08、N09、all）。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `simulation`: 一時出力先へ保存場を作る関数。
- `analysis`: 保存場を読み、解析結果を作る関数。
- `plotting`: 保存場と解析結果から図を作る関数。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

nothing。選択IDの各6ファイルを反映し、他IDと発展出力を保持する。生成失敗は既存出力を保持し、反映失敗は復元する。
"""
function main(;
    selection = "all",
    output_dir = N08N09Simulation.DEFAULT_OUTPUT_DIR,
    simulation = N08N09Simulation.main,
    analysis = N08N09Analysis.main,
    plotting = N08N09Plots.main,
    publish_options...,
)
    N08N09Simulation.selection_ids(selection)

    # 一時領域で処理を完了してから、公式出力を反映する。
    mktempdir() do stage
        simulation(; selection, output_dir = stage)
        analysis(; selection, input_dir = stage, output_dir = stage)
        plotting(; selection, input_dir = stage, output_dir = stage)
        N08N09Simulation.check_complete(stage; selection)
        N08N09Simulation.publish(
            stage,
            output_dir,
            N08N09Simulation.output_names(selection);
            publish_options...,
        )
    end
    println("$selection の公式出力を検証してまとめて反映しました。")
end
abspath(PROGRAM_FILE) == (@__FILE__) &&
    main(; selection = N08N09Simulation.cli_selection(ARGS))
end
