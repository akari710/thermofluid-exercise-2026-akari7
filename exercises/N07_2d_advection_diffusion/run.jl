# N07の提供済み一括入口。数値・解析の編集先は次のファイルです。
# src/N07Transport.jl: thermal_stable_timestep / thermal_fluxes! / thermal_step!と、burgers_stable_timestep / burgers_step!の5関数。
# analyze.jl: heat_budget。tests.jl: 配布必須テストを読み、「自作」枠に自分のテストを書きます。
# 読む順: 課題N07（面添字・熱収支） → 数値5関数とtests.jl → simulate.jl → analyze.jl → plot.jl → このmain。
# 温度・Burgersの保存場と境界輸送を計算し、保存データを解析・作図します。
# 計算 → 解析 → 作図 → 完了検査を一時領域で終え、公式8出力へ反映します。

module N07Run
include("simulate.jl")
include("analyze.jl")
include("plot.jl")

"""
    main(;
        output_dir = N07Simulation.DEFAULT_OUTPUT_DIR,
        simulation = N07Simulation.main,
        analysis = N07Analysis.main,
        plotting = N07Plots.main,
        publish_options...,
    )

計算・解析・作図・完了検査を一時領域で終え、全成果物を反映する。
既定実行には `src/N07Transport.jl` の数値5関数と `analyze.jl` の `heat_budget` の実装が必要。
数値・解析TODOの未実装エラーでは反映へ進まず、既存出力を保持する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `simulation`: 一時出力先へ保存場を作る関数。
- `analysis`: 保存場を読み、解析結果を作る関数。
- `plotting`: 保存場と解析結果から図を作る関数。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

nothing。公式8出力を反映する。どの生成段階の失敗でも既存出力を保持し、反映失敗は復元する。
"""
function main(;
    output_dir = N07Simulation.DEFAULT_OUTPUT_DIR,
    simulation = N07Simulation.main,
    analysis = N07Analysis.main,
    plotting = N07Plots.main,
    publish_options...,
)
    mktempdir() do stage
        simulation(; output_dir = stage)
        analysis(; input_dir = stage, output_dir = stage)
        plotting(; input_dir = stage, output_dir = stage)
        N07Simulation.check_complete(stage)
        N07Simulation.publish(
            stage,
            output_dir,
            N07Simulation.OUTPUT_NAMES;
            publish_options...,
        )
    end
    println("N07の公式8出力をまとめて反映しました。")
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
