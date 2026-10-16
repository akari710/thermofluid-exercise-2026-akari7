# 時間積分の比較: 今年度は任意。実装するのはrk2_stepのTODO 1箇所。
# 読む順: このファイル → provided_methods.jl → tests.jl → provided_support.jl。
# 既存N03の実装や進捗とは独立し、実験・描画・保存は配布済み。
module TimeIntegrationExample

using LinearAlgebra, TOML

# 先に配布手法と入力検証を読み込み、その後に学生のRK2を定義する。
include("provided_methods.jl")

"""
    rk2_step(f, y, t, dt)

陽的中点法で、状態を時刻 `t` から `t + dt` へ1ステップ進める。

# 引数

- `f`: 状態と時刻を受け取る右辺関数 `f(y, t)`。
- `y`: 有限な実数スカラー（`Bool` を除く）、または空でない有限な `Vector{Float64}`。
- `t`: 現在の有限な時刻。
- `dt`: 有限な正の時間刻み。

不正な状態・時刻・時間刻みは、計算前に `ArgumentError` で拒否する。
配布済みの `rhs(f, y, t)` は入力状態を保護し、右辺の形と有限性を確認する。

# 返り値

実装後は、次時刻の状態を返す。
入力 `y` を変更せず、配列の場合は `y` と独立した配列を返す。

# 受講生のToDo

区間の中点での傾きを使う、陽的中点法の二段更新を実装する。
配布時は `未実装 N03: rk2_step (任意比較)` のエラーで停止する。
"""
function rk2_step(f,y,t,dt)
    # 未実装の状態でも、不正な入力は更新前に検出する。
    validate_step(y,t,dt)

    # TODO(N03): 陽的中点法の二段更新を実装する（今年度は任意比較）。
    error("未実装 N03: rk2_step (任意比較)")
end

# 定義済みの手法を使う実験・診断・描画・保存を最後に読み込む。
include("provided_support.jl")
end

# 直接実行したときだけmainを呼ぶ。includeでは関数の定義だけを行う。
if abspath(PROGRAM_FILE) == @__FILE__
    # --cnは任意比較の指定。それ以外の引数は出力先1つまでとする。
    arguments = filter(!=("--cn"),ARGS)
    length(arguments)<=1 || throw(ArgumentError("使い方: run.jl [出力先] [--cn]"))

    # 指定がなければ学生リポジトリ内のscratchへ保存する。
    directory = isempty(arguments) ? TimeIntegrationExample.default_output_dir() : only(arguments)
    TimeIntegrationExample.main(;output_dir=directory,include_cn="--cn" in ARGS)
end
