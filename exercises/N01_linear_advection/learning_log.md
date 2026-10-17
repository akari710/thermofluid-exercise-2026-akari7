# N01 学習ログ

## 予想

- 風上差分と中心差分＋陽Eulerの波形：風上差分は，時間の経過とともに+x方向に移動し，矩形がなくなってくる．中心差分+陽Eulerは，矩形の形を保ったまま移動する．
- 安定性を判断する数値：courant

## 式とコード

- 矩形初期条件とコードの対応：
function rectangular_initial_condition(
    x::AbstractVector{<:Real};
    base::Real=1.0,
    plateau::Real=2.0,
    plateau_start::Real=0.5,
    plateau_end::Real=1.0,
)
    validate_initial_condition_inputs(x, base, plateau, plateau_start, plateau_end)

    #=
    TODO(N01):
    矩形状の初期分布を実装する。
    各座標xiがplateau_start <= xi <= plateau_endを満たすか判定する。
    =#
    u = similar(x)
    for i in eachindex(x)
        xi = x[i]
        if plateau_start <= xi <= plateau_end
            u[i] = plateau
        else
            u[i] = base
        end
    end
    return u
end
- 風上差分・中心差分の式と添字の対応：添え字iは1から始まる．有効添字は，後退差分が i = 2, ..., length(u)，中心差分はがi = 2, ..., length(u)-1である．
風上差分
u_new[i] = u_old[i] - courant * (u_old[i] - u_old[i-1])
中心差分
u_new[i] = u_old[i] - courant * (u_old[i+1] - u_old[i-1]) / 2
- 境界条件、時間ループ、バッファ交換の説明：【境界条件】正の速度では情報が左から領域へ入り，右へ出る．この1次元線形移流課題では左端を1に固定し，右端は左隣の値をコピーしてゼロ勾配を表す．関数内では，u[1] = 1.0，u[end] = u[end - 1]のように書く．1ステップの内部点の更新後に，この境界条件を適用する．【時間ループ】局所的な1ステップの更新を繰り返すと，時間発展を計算できる．時間発展ループは
for _ in 1:steps
    step!(u_new, u_old, c, dt, dx)
    apply_boundary!(u_new)
    u_old, u_new = u_new, u_old
end
のように書かれ，各段階では，古い配列だけを参照して新しい配列の内部点を計算し，境界条件を適用してから次の段階へ進む．【バッファ交換】u_old, u_new = u_new, u_oldで二つの変数名の参照先を交換する．計算済みのu_newを次の段階のu_oldとして使っている．

## 結果と解釈

- `nx`、`dx`、`steps`、`dt`、実効CFL：nx=81，dx=0.025，steps=40，dt=0.0125，実効CFL=0.5
- 二手法の最小・最大値とオーバーシュート・アンダーシュート：【風上差分】最小1.0，最大1.9993204517450067，オーバーシュート0.0，アンダーシュート0.0　【中心差分】最小-10.420010468371677，最大13.492726642559711，オーバーシュート11.492726642559711，アンダーシュート11.420010468371677
- 図：[風上差分](results/upwind.png)、[中心差分＋陽Euler](results/centered-euler.png)
- 数値記録：[summary.toml](results/summary.toml)
- 予想との比較、解釈、限界：風上差分は概ね予想通りであった．中心差分はオーバーシュート，アンダーシュートいずれも大きな値を取り，予想と違って負の値も取っていた．

## 自分のテスト

- 入力と手計算した期待値：自作テストは以下のようである．
    # 1. 矩形状初期分布
    x = [0.0, 0.5, 1.0, 1.5, 2.0]

    u = rectangular_initial_condition(
        x;
        base=1.0,
        plateau=3.0,
        plateau_start=0.5,
        plateau_end=1.5
    )

    @test u == [1.0, 3.0, 3.0, 3.0, 1.0]
    # 2. 風上差分法
    u_old = [1.0, 2.0, 4.0, 3.0, 1.0]
    u_new = similar(u_old)

    upwind_step!(u_new, u_old, 1.0, 0.5, 1.0)

    @test u_new == [1.0, 1.5, 3.0, 3.5, 1.0]
    # 3. 中心差分法
    u_new = similar(u_old)

    centered_step!(u_new, u_old, 1.0, 0.5, 1.0)

    @test u_new == [1.0, 1.25, 3.75, 3.75, 1.0]
    # 4. 境界条件
    u = [10.0, 2.0, 3.0, 4.0, 20.0]

    apply_boundary!(u; left_value=1.0)

    @test u == [1.0, 2.0, 3.0, 4.0, 4.0]
    # 5. 安定性の比較
    upwind = simulate(
        scheme=:upwind,
        nx=81,
        c=1.0,
        cfl=0.5,
        t_final=0.5
    )

    centered = simulate(
        scheme=:centered,
        nx=81,
        c=1.0,
        cfl=0.5,
        t_final=0.5
    )

    @test minimum(upwind.u) >= 1.0 - 1e-10 &&
        maximum(upwind.u) <= 2.0 + 1e-10 &&
        (minimum(centered.u) < 1.0 - 1e-6 ||
            maximum(centered.u) > 2.0 + 1e-6)
- 保証すること・保証しないこと：
- 実行コマンドと結果：【julia --project=. exercises/N01_linear_advection/run.jl】N01の出力を書き込みました: /home/akari710/thermofluid-exercise-2026-akari7/exercises/N01_linear_advection/results，【julia --project=. -e 'using Pkg; Pkg.test()'】...Testing ThermofluidExercise tests passed
- 必須TODOと自作テストの入力・期待値を選んだ理由：1桁の整数の方が計算しやすいため．

## AI利用・出典

- 依頼内容（利用なしの場合は「利用なし」）：テストを実行したところ，
N01 自作テスト: Error During Test at /home/akari710/thermofluid-exercise-2026-akari7/exercises/N01_linear_advection/tests.jl:75
  Got exception outside of a @test
  UndefVarError: `rectangular_initial_condition` not defined in `Main`
  Suggestion: check for spelling errors or missing imports.
  と出てきてしまい，どこを修正すればよいのか教えてほしい
- 重要な提案：test.jlの冒頭にusing .N01LinearAdvectionを追加するとよい．
- 採用・修正・却下と理由：run.jlでは，rectangular_initial_conditionがN01LinearAdvectionというモジュールに属しているが，tests.jlではrectangular_initial_condition(...)と直接呼び出しており，Juliaが関数を見つけられなくなっていた．
- 採用素材の出典、変更点、利用条件（なければ「なし」）：なし


## 理解度チェック・LETUS提出

- 対応する授業ID：N01
- LETUS提出日：
- 提出済み確認：
- 理解できた点：
- 残った疑問：
- 対話全文はLETUSへ提出し，このリポジトリには含めていない：
