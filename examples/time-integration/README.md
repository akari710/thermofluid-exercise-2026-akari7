# 時間積分の比較（2026年度・任意）

B3向けに，常微分方程式で時間次数を確認し，同じ拡散問題で安定性・精度・計算量を比較する独立例です．
既存N03の基準は前進Euler＋空間2次中心差分です．
今年度のN03完了条件と通常の`Pkg.test()`にはこの例を含めません．
次年度のN03必修化は別作業です．

学生が編集するのは`run.jl`の`rk2_step(f,y,t,dt)`のTODO 1箇所です．
陽的中点法の二段更新を実装してください．
前進Euler・古典的RK4・後退Euler・任意のCrank–Nicolsonは`provided_methods.jl`の配布実装を読み，式と結果を比較してください．
後退Eulerは旧値から新値を求める線形方程式を解き，固定dtでは行列分解を再利用します．
行列組立てを学生の必須実装にはしません．
既存N03や共通パッケージをinclude/importする必要はありません．

## APIと実行

右辺は`f(y,t)`です．
状態は有限実数スカラーまたは空でない`Vector{Float64}`，時刻は有限，dtは有限正とします．
入力状態を変更せず，配列の返り値は入力と独立させてください．
入力検証と右辺呼出しの検証`rhs(f,y,t)`は配布済みです．
未実装時にはテストと実行が`未実装 N03: rk2_step (任意比較)`で失敗します．
仮の結果で通さず，TODOを完成してから次のコマンドを学生リポジトリrootで実行してください．
既存のJulia環境を使います．

```fish
julia --project=. examples/time-integration/tests.jl
julia --project=. examples/time-integration/run.jl
```

既定出力は`scratch/time-integration/`の`time-convergence.png`，`diffusion-comparison.png`，`summary.toml`です．
今年度の公式提出物には含めません．
出力先の指定と任意CNの追加は次のように実行します．

```fish
julia --project=. examples/time-integration/run.jl scratch/time-integration-cn --cn
```

includeだけでは計算・表示・保存・進捗更新を行いません．
入力不正は状態・出力を書き換える前に`ArgumentError`とします．
未実装・計算失敗・描画失敗・TOML保存失敗では旧出力を保持します．
一時領域で3ファイルを完成させてから反映し，反映途中の失敗では旧3ファイルに戻します．
無関係なファイルと`course_progress.toml`は変更しません．

## 比較する結果

実験Aは`y'=-y, y(0)=1, T=1`を`N=8,16,32,64`，`dt=T/N`で積分します．
`abs(y_N-exp(-1))`と隣接誤差比の`log2`を読み，前進Euler・RK2・RK4・後退Eulerの観測次数1・2・4・1を比べてください．
ここには空間誤差がありません．
積分中のRHS評価回数，分解回数，線形求解回数を比べます．
診断用計算をこれらの回数に含めず，実時間を採点基準にしません．

実験Bは`[0,2]`，`D=0.1`，両端を含む41点，`dx=0.05`，両端Dirichlet 0の空間2次中心差分です．
内部39点を未知数とし，`phi_k(j)=sin(k*pi*j/40)`，`j=0,…,40`として初期値を`phi_1+0.01*phi_39`とします．
`T=0.1`に対し`N=128`（Fo=0.03125）と`N=4`（Fo=1）を比べます．
各Nで全手法の格子・初期値・dtを揃え，陽解法のためにdtを縮小しません．
参照は`lambda_k=-4D/dx^2*sin(k*pi/80)^2`による半離散厳密解`exp(lambda_1*t)*phi_1+0.01exp(lambda_39*t)*phi_39`です．
同じ空間離散化に対する時間誤差を見るため，連続PDEの厳密解を参照にしません．
最終分布，最大絶対誤差，高周波モード39の係数履歴を比較してください．
粗いdtで陽解法の高周波が有限のまま増幅する結果も保存します．
非有限値は手法・ステップ・時刻・dtを示して停止します．
任意CNでは高周波係数の符号交代も見られます．
この線形拡散問題での無条件安定は，任意のdtで正確または非振動であることを意味しません．

`summary.toml`の`schema_version=1`，`example_id="time-integration"`がこの例を識別します．
`time_convergence`には手法ごとのN・dt・実終端時刻・最終値・誤差・全3区間の次数・計算量を記録します．
`diffusion_problem`には格子・D・境界・初期モードと係数・固有値・参照の種類を，`diffusion`にはNと手法ごとのFo・時刻・モード係数履歴・初期／最終／参照分布・誤差・計算量を記録します．

授業と既存課題は[公開N03授業](https://t2lab-it.github.io/thermofluid-exercise-2026/lessons/N03.html)と[公開N03課題](https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/N03.html)を参照してください．
時間積分の公開解説・導線と，必要なexamplesだけの取得・復元手順はセッション2で整備します．
