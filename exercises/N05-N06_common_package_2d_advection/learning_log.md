# N05-N06 学習ログ

## N05：抽出前の予想とbaseline

- 残す挙動・実際の重複・N04選択：共通化しても，これまでと同じ引数で関数を呼び出せて，同じ条件なら同じ計算結果が得られるようにしたいと考えた．特に，周期境界で端の隣をどの点にするかと，更新中に旧配列を書き換えないことに気を付けた．N02とN04では周期境界の隣の点を求める処理があり，N01〜N04には新旧配列のチェックや，最終時刻に合わせて刻みを調整する処理があった．N04では線形モデル（linear）を選択した．
- baselineを取得したcommitと対象コード：共通化する前に N05.jl baseline を実行して，その時点のN01〜N04の計算結果を baseline.toml に保存した．基準を取得したときのコミットIDは"a0b60ccd033696f135042ed782ee341affd9ba33"である．

## N05：採用したdiffと回帰

- 共通APIへ移した処理と、課題に残した処理：周期境界で左右の隣を求める処理，線形流束とBurgers流束，新旧配列のチェック，時間刻みの調整を共通パッケージにまとめた．各課題の更新式や境界処理は，共通パッケージ内の対応するモジュールに移した．一方，初期条件，計算条件，結果の集計や作図は各課題に残した．N04で選んだモデルも課題側で設定し，共通側の関数に引数として渡すようにした．
- テスト・代表結果の差とその理由：N05.jl verify を実行し，共通化の前後でN01〜N04の代表結果を比較した．比較結果は一致した．regression.tomlでは，基準ファイルと現在のコードのSHA-256，実行時のコミットID，N04の選択モデルを照合した．過去課題のテスト結果はすべて合格だった．
- 復帰点を使った場合の採用差分と説明：復帰点は使用していない．

## N06：二次元更新と数値検証

- 軸・周期端・角・方向別寄与・合成CFL：u_old[i, j] の i をx方向，jをy方向に対応させた．各点を更新するときはx方向とy方向の流束の差を，どちらも更新前の u_old から求めた．周期境界では先頭の一つ手前を末尾とするので，角の (1, 1) では，x方向に (nx, 1)，y方向に (1, ny) を参照する．時間刻みは，二方向のCFL数の合計C_x+C_y<=1となるように決めた．
- 複数ステップの保存量、自作テスト：値が場所によって異なる非対称な3×4の配列を使い，dx = 0.5，dy = 0.25，dt = 0.1，cx = 1.0，cy = 0.75の条件で10ステップ進めた．初期積分は，更新関数を使わずに，全要素の和78にdx　*　dyを掛けて9.75と求めた．毎ステップの積分を9.75と比べ，差が9.75*10^{-12}以下かを確認した．許容誤差は，繰り返し計算で生じる丸め誤差を考えて，初期積分の大きさに10^{-12}を掛けて決めた．テストは成功した．
- 3格子の誤差・観測次数と考察：公式の3格子でt=1まで計算し，同じ時刻の解析解と比べた．各格子点で数値解と解析解の差を二乗し，その平均の平方根をL_2誤差とした．40×30，80×60，160×120の順に，誤差は0.04586691138478507, 0.02406462411492619, 0.012332831539922874だった．格子幅を半分にしたときの誤差の比から求めた観測次数は，順に0.9305398727322467, 0.9644098174384694だった．

## N06：保存場の再解析と再作図

- 分散を追加した結果、HDF5の変更前後のSHA-256：analyze.jl に，保存された各時刻について空間分散を求める処理を追加した．まず全格子点の平均を求め，次に各点の平均との差を二乗して足し，最後に全格子点数で割った．summary.tomlの diagnostics_completeは，分散を追加する前がfalse，追加した後がtrueだった．
- 表示設定の変更、入力HDF5・summaryのSHA-256：plot.jlでは，授業ページの画像と同様の画像が得られたことから特に変更しなかった．plots.tomlに記録された入力ファイルのSHA-256は，fields.h5が"3c9350f39d725d55331c140b5d269e9dad613745f681701485dde3e444a261bc"，summary.tomlが"be40a4d61b4c5a48055a094c5ab818c46e80a96c0fe89d0ec400f1d4390dc9d4"だった．
- 保存されなかった時刻からは求められない量：保存した時刻については，その時刻の分布や空間分散を計算できる．例えばt=0.25とt=0.50の間で，数値解がいつ最大になったかを正確に知るには，その間の計算結果も保存しておく必要がある．

## 理解度チェック・提出前確認

- N05・N06それぞれの理解度チェックと振り返り：理解度チェックは提出済み．
- 実行コマンド、テスト結果、成果物、未解決事項：実行コマンドは「git rev-parse HEAD」「julia --project=. exercises/N05-N06_common_package_2d_advection/N05.jl baseline」「julia --project=. exercises/N05-N06_common_package_2d_advection/N05.jl verify」「julia --project=. exercises/N05-N06_common_package_2d_advection/simulate.jl」「julia --project=. exercises/N05-N06_common_package_2d_advection/analyze.jl」「julia --project=. exercises/N05-N06_common_package_2d_advection/plot.jl」「julia --project=. -e 'using SHA; println(bytes2hex(sha256(read("exercises/N05-N06_common_package_2d_advection/results/N06/fields.h5"))))'」「julia --project=. exercises/N05-N06_common_package_2d_advection/run.jl」「julia --project=. exercises/N05-N06_common_package_2d_advection/tests.jl」「julia --project=. -e 'using Pkg; Pkg.test()'」である．テストは成功．N05の baseline.toml と regression.toml，N06の fields.h5，summary.toml，図3枚，plots.tomlは"3c9350f39d725d55331c140b5d269e9dad613745f681701485dde3e444a261bc"である．
