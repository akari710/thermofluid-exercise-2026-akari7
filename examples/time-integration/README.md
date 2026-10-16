# 時間積分の比較（2026年度・任意）

B3向けに，常微分方程式で時間精度の次数を求め，拡散問題で終時刻の誤差・高周波成分の増減・処理回数を比べる例です．
既存N03の基準は前進Euler＋空間2次中心差分です．

学生が編集するのは`run.jl`の`rk2_step(f,y,t,dt)`のTODO 1箇所です．
陽的中点法の二段更新を実装してください．
前進Euler・古典的RK4・後退Euler・任意のCrank–Nicolsonは`provided_methods.jl`の配布実装を読み，式と結果を比較してください．
後退Eulerは新時刻の状態を未知数とする連立一次方程式を解き，手法と時間刻みを決めるごとに係数行列を1回LU分解して全ステップで再利用します．

## 新規フォルダだけの取得

公開解説は[時間積分の精度・安定性・計算量](https://t2lab-it.github.io/thermofluid-exercise-2026/advanced/numerical-schemes.html)を参照してください．
既存の`upstream`が授業の学生配布元を指すことを，共通の[配布元の更新手順](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/workflow.html#material-updates)で確認します．
学生リポジトリの最上位ディレクトリの，現在の課題用ブランチで次のブロック全体を実行します．

```bash #time-integration-fetch
bash -c '
target=examples/time-integration
root=$(git rev-parse --show-toplevel) || exit 1
if [ "$(pwd -P)" != "$root" ]; then
    echo "学生リポジトリのrootで実行してください．" >&2
    exit 1
fi
if [ -e "$target" ] || [ -L "$target" ] ||
   [ -n "$(git ls-files -- "$target")" ]; then
    echo "同名パスの内容と差分を確認してください．" >&2
    exit 1
fi
if [ -L examples ] || { [ -e examples ] && [ ! -d examples ]; }; then
    echo "examplesの状態を確認してください．" >&2
    exit 1
fi
git fetch --no-tags upstream main || exit 1
if [ "$(git cat-file -t "FETCH_HEAD:$target")" != tree ]; then
    echo "配布元のexamples/time-integrationを確認してください．" >&2
    exit 1
fi
git restore --source=FETCH_HEAD --worktree -- "$target"
'
```

これは初回の取得手順です．
既存の課題とステージ領域を保ち，新規フォルダだけが未追跡で加わります．
`git status --short`と`git diff --cached`で取得結果と元の作業を確認します．
通常の配布元更新は[共通の更新手順](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/workflow.html#material-updates)を参照してください．

## APIと実行

右辺関数`f(y,t)`は状態`y`と時刻`t`から時間微分を返し，`rk2_step(f,y,t,dt)`は時間刻み`dt`だけ進めた新しい状態を返します．
状態は有限実数のスカラーまたはベクトル（1要素以上の`Vector{Float64}`），時刻は有限実数，時間刻み`dt`は有限正とします．
入力状態を保ち，配列の返り値は入力から独立させてください．
入力検証と右辺呼出しの検証`rhs(f,y,t)`は配布済みです．
実装後は，既存のJulia環境を使い，学生リポジトリの最上位ディレクトリで次の順に実行します．

```bash
julia --project=. examples/time-integration/tests.jl
julia --project=. examples/time-integration/run.jl
```

既定出力は`scratch/time-integration/`の`time-convergence.png`，`diffusion-comparison.png`，`summary.toml`です．
出力先の指定と任意CNの追加は次のように実行します．

```bash
julia --project=. examples/time-integration/run.jl scratch/time-integration-cn --cn
```

## 比較する結果

実験Aでは方程式`y'=-y`，初期値`y(0)=1`，最終時刻`T=1`を固定し，各手法で`N=8,16,32,64`，`dt=T/N`として時間刻みを半分ずつにします．
終時刻の絶対誤差`E_N=abs(y_N-exp(-1))`と，2計算の誤差から求める観測次数`log2(E_N/E_2N)`を読みます．
前進Euler・RK2・RK4・後退Eulerの順に，観測次数が1・2・4・1へ近づくことを確認してください．
同じNでは時間積分法だけを変え，絶対誤差と処理回数を比べます．

処理回数として数えるのは，1回の積分の準備と全Nステップで行う右辺関数の評価（RHS評価）・行列分解・連立一次方程式の求解です．
`cost`の`rhs_evaluations`・`factorizations`・`linear_solves`は，それぞれの関数が呼ばれた実際の回数です．
後退Eulerでは1回のLU分解とN回の求解を数えます．
右辺評価と線形求解の負担は行列の大きさや構造に応じて異なるため，3種類の回数を個別に比べます．

実験Bは`[0,2]`，拡散係数`D=0.1`，両端を含む41点，格子幅`dx=0.05`の空間2次中心差分で，両端の値を0に固定します．
内部39点を未知数とし，`phi_k(j)=sin(k*pi*j/40)`，`j=0,…,40`として初期値を`phi_1+0.01*phi_39`とします．
最終時刻`T=0.1`に対し`N=128`（`dt=0.00078125`，Fo=0.03125）と`N=4`（`dt=0.025`，Fo=1）を比べます．
同じNでは格子・拡散係数・境界値・初期値・最終時刻・時間刻みを揃え，時間積分法だけを変えます．
刻みの影響は，各手法でN=128とN=4を比べて調べます．
参照は`lambda_k=-4D/dx^2*sin(k*pi/80)^2`による半離散厳密解`exp(lambda_1*t)*phi_1+0.01*exp(lambda_39*t)*phi_39`です．
同じ半離散系を基準に，時間積分だけの誤差を比較します．
T=0.1での分布と，内部39点での最大絶対誤差`maximum(abs.(u - u_ref))`を比較してください．
高周波モード39は隣り合う格子点で符号が交代する分布で，その係数を内部39点の`sum(u .* phi39) / sum(phi39 .^ 2)`で求めます．
初期時刻と各ステップ後の時刻`n*dt`（n=0,…,N）で，係数を半離散厳密解の`0.01*exp(lambda_39*n*dt)`と比べます．
N=4での陽解法の係数の増幅と，任意のCrank–Nicolsonでの符号交代を確認します．
後退EulerではN=4での係数の減衰を確認し，N=4とN=128の最終誤差も比べます．
増幅率と誤差の関係は[公開解説](https://t2lab-it.github.io/thermofluid-exercise-2026/advanced/numerical-schemes.html#backward-euler)を参照してください．

`summary.toml`の`time_convergence`でNごとの絶対誤差・処理回数と，N=8→16・16→32・32→64の観測次数を確認できます．
`diffusion`の`max_abs_error`は最終分布の最大絶対誤差，`high_mode_coefficients`はモード39の係数を時刻順に並べた配列です．

授業と既存課題は[公開N03授業](https://t2lab-it.github.io/thermofluid-exercise-2026/lessons/N03.html)と[公開N03課題](https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/N03.html)を参照してください．
