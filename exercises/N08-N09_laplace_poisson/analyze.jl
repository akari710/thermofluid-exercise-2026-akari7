# N08/N09: 保存場の読取り → 独立残差・解析解誤差・格子収束 → 診断保存。

module N08N09Analysis
include("provided_support.jl")

"""
    main(;
        selection = "all",
        input_dir = DEFAULT_OUTPUT_DIR,
        output_dir = DEFAULT_OUTPUT_DIR,
        publish_options...,
    )

選択した保存場から独立残差・誤差・収束を診断する。

# 引数

- `selection`: 処理する内容ID（N08、N09、all）。
- `input_dir`: 保存場HDF5を読むディレクトリ。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

nothing。run_id混在や診断失敗は拒否し、選択IDのsummary.tomlを復元付き反映する。
"""
function main(;
    selection = "all",
    input_dir = DEFAULT_OUTPUT_DIR,
    output_dir = DEFAULT_OUTPUT_DIR,
    publish_options...,
)
    ids = selection_ids(selection)
    pairs = Dict(
        id => read_fields(joinpath(input_dir, id, "fields.h5"); expected_id = id) for
        id in ids
    )

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        length(unique(pairs[id].metadata["run_id"] for id in ids)) == 1,
        "run_idが混在しています",
    )

    # 一時領域で処理を完了してから、公式出力を反映する。
    staged(
        output_dir,
        Tuple(joinpath(id, "summary.toml") for id in ids);
        publish_options...,
    ) do stage
        for id in ids
            mkpath(joinpath(stage, id))
            write_toml(
                joinpath(stage, id, "summary.toml"),
                diagnostics(pairs[id], file_sha(joinpath(input_dir, id, "fields.h5"))),
            )
        end
    end
    println("$(join(ids,"・")) 独立残差・解析解誤差・3格子収束を解析しました。")
end
abspath(PROGRAM_FILE) == (@__FILE__) && main(; selection = cli_selection(ARGS))
end
