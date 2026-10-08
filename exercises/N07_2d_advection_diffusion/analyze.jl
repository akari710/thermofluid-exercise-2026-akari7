# N07: 保存場読取り → 受講生の熱収支TODO → 診断保存。

module N07Analysis
include("provided_support.jl")

"""
    heat_budget(time, heat, advective_integrals, diffusive_integrals)

流入を正とする辺別区間積分から、熱量変化と収支残差を再構成する。

# 引数

- `time`: 0から始まる有限な狭義単調増加の保存時刻列。長さntは2以上。
- `heat`: 時刻に対応する有限な熱量列。
- `advective_integrals`: 流入を正とする有限な移流熱輸送の(4,nt-1)行列。行順はwest,east,south,north、列は隣り合う保存時刻の区間。
- `diffusive_integrals`: 同じ辺順・区間順の有限な拡散熱輸送の(4,nt-1)行列。

入力4配列は変更しない。不正入力は提供のvalidate_budgetでArgumentError。

# 返り値

実装後は `net_input(nt-1)`、`cumulative_input(nt)`、`residual(nt)` のNamedTuple。
後二つの先頭は0。residualは初期熱量からの変化から累積流入を引いた列。

# 受講生のToDo

保存区間ごとに4辺の移流・拡散積分を合計し、正味流入とその累積、初期熱量からの変化との差を再構成する。
入力検証は提供済み。積分はsimulate.jlが各stepの辺別レート×刻みを保存区間に累積した値なので、ここで刻みを再び掛けない。
[N07課題「辺別熱輸送と熱収支」](https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/N07.html#heat-budget)を参照する。配布状態では検証後に未実装エラーで停止する。
"""
function heat_budget(time, heat, advective_integrals, diffusive_integrals)
    validate_budget(time, heat, advective_integrals, diffusive_integrals)
    # TODO(N07): 各区間の正味流入、累積流入、初期熱量からの変化との差を求める。
    error("未実装 N07: heat_budget")
end

"""
    main(;
        input_dir = DEFAULT_OUTPUT_DIR,
        output_dir = DEFAULT_OUTPUT_DIR,
        budget = heat_budget,
        publish_options...,
    )

同じ実行の両保存場を読み、熱収支と数値診断を保存する。

# 引数

- `input_dir`: 保存場HDF5を読むディレクトリ。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `budget`: 熱量と辺別区間積分から熱収支を返す学生関数。
- `publish_options`: 容量検査・反映・復元へ渡すキーワード引数。

# 返り値

解析辞書。source_sha256を記録してsummary.tomlを復元付きで反映する。熱収支TODOの失敗では既存解析を保持する。
"""
function main(;
    input_dir = DEFAULT_OUTPUT_DIR,
    output_dir = DEFAULT_OUTPUT_DIR,
    budget = heat_budget,
    publish_options...,
)
    # 保存済みの値と出自を検証してから、解析や作図に使う。
    data = read_pair(input_dir)

    # 計算値と条件を、保存する診断情報にまとめる。
    doc = diagnostics(data, budget)
    doc["source_sha256"] = input_hashes(input_dir)

    # 一時領域で処理を完了してから、公式出力を反映する。
    staged(output_dir, ("summary.toml",); publish_options...) do stage
        write_toml(joinpath(stage, "summary.toml"), doc)
    end
    println("N07の保存場と区間熱輸送を解析しました。plot.jlを再実行してください。")
    doc
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
