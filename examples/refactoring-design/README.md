# 共通化の設計判断：コード一式

[発展資料「共通化の設計判断」](https://t2lab-it.github.io/thermofluid-exercise-2026/advanced/refactoring-design.html)に掲載したコードです．
二つの境界処理と入力検証を共通モジュールへ移し，共通化の前後と比較案を実行できます．

## ファイル構成

| ファイル | 役割 |
|---|---|
| `before.jl` | 共通化前の境界処理と入力検証 |
| `common.jl` | 境界処理と入力検証をまとめた共通モジュール |
| `after.jl` | 共通化後の呼出し用関数 |
| `wide.jl` | 境界の種類を一つの関数で選ぶ比較案 |
| `strict.jl` | 必要な配列長を変えた比較案 |
| `checks.jl` | 共通化前後を同じ条件で調べる検査 |
| `compare.jl` | 比較案の動作を調べる検査 |

## 実行方法

学生リポジトリのルートから移動します．

```sh
cd examples/refactoring-design
```

次の各コマンドを実行します．
`checks.jl`は抽出前・抽出後・共通側をそれぞれ検査し，`compare.jl`は二つの比較案を検査します．

```sh
julia --startup-file=no --compiled-modules=no --project=. checks.jl before
julia --startup-file=no --compiled-modules=no --project=. checks.jl after
julia --startup-file=no --compiled-modules=no --project=. checks.jl common
julia --startup-file=no --compiled-modules=no --project=. compare.jl
```

それぞれ新しいJuliaプロセスで実行されます．
期待した結果を確認できると，全検査が成功し，終了コードは0になります．
