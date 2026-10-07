# N05/N06: 回帰確認 → simulate → analyze → plot → 全出力の反映を一時領域で行います。

module N0506Run
include("N05.jl")
include("simulate.jl")
include("analyze.jl")
include("plot.jl")

"""
    main(;
        output_dir = joinpath(@__DIR__, "results"),
        baseline_path = joinpath(output_dir, "N05", "baseline.toml"),
        simulation = N06Simulation.main,
        analysis = N06Analysis.main,
        plotting = N06Plots.main,
        publish_options...,
    )

数値回帰・計算・解析・作図を一時領域で順に完了して全出力を反映する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `baseline_path`: 抽出前の数値結果を保存したbaseline TOMLのパス。
- `simulation`: 一時出力先へ保存場を作る関数。
- `analysis`: 保存場を読み、解析結果を作る関数。
- `plotting`: 保存場と解析結果から図を作る関数。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

nothing。N05/regression.tomlとN06公式6ファイルを反映し、baselineは保持する。どの段階の失敗でも既存出力を保ち、反映失敗は復元する。
"""
function main(;
    output_dir = joinpath(@__DIR__, "results"),
    baseline_path = joinpath(output_dir, "N05", "baseline.toml"),
    simulation = N06Simulation.main,
    analysis = N06Analysis.main,
    plotting = N06Plots.main,
    publish_options...,
)
    mktempdir() do stage
        N05Regression.verify(; baseline_path, output_dir = joinpath(stage, "N05"))
        target = joinpath(stage, "N06")
        simulation(; output_dir = target)
        analysis(; input_path = joinpath(target, "fields.h5"), output_dir = target)
        plotting(;
            input_path = joinpath(target, "fields.h5"),
            summary_path = joinpath(target, "summary.toml"),
            output_dir = target,
        )
        names = vcat(
            [joinpath("N05", "regression.toml")],
            [joinpath("N06", n) for n in N06Simulation.OUTPUT_NAMES],
        )
        N06Simulation.publish(stage, output_dir, names; publish_options...)
    end
    println("N05回帰とN06の全出力を反映しました。baselineは保持しています。")
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
