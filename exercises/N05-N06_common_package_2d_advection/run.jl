# N05/N06の提供済み一括入口。編集先はこのmainではなく、次の数値・解析関数です。
# N05: src/ThermofluidExercise.jlの共通APIとN01〜N04の移植先を実装し、各run.jlの対応関数を共通APIへ委譲します。
# N05.jlは移植前のbaseline取得と移植後の比較の提供入口です。baselineを保持して比較します。
# N06: src/N06Advection.jlのstable_timestep / advection_step!、analyze.jlのspatial_varianceを実装します。
# tests.jlの配布必須テストを読み、同じファイルの「自作」枠に自分のテストを書きます。
# 読む順: 課題N05/N06 → 上記の関数とtests.jl → simulate.jl → analyze.jl → plot.jl → このmain。
# 保存の流れ: simulate.jl → fields.h5 → analyze.jl → summary.toml → plot.jl → 図・plots.toml。
# 解析と作図は保存場を読み、数値場を再計算しません。一括入口は全段階を一時領域で完了してから反映します。

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
既定実行にはN05の共通化・旧入口の委譲とbaseline、N06の数値2関数と `spatial_variance` の実装が必要。
数値TODOへ到達すると未実装エラーで停止する。
分散だけ未完成なら `analyze.jl` の単独実行は基本診断を保存し、`diagnostics_complete=false` とする。
この一括入口はその後の作図で失敗するため、公式出力の反映へ進まない。

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
