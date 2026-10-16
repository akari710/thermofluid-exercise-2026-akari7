# 入力検証と時間積分手法の配布実装。学生の必須編集箇所ではありません。
# run.jlから読み込み、TimeIntegrationExampleモジュールの一部になります。
# 読む順: 入力検証 → 陽解法 → 線形問題の陰解法 → 拡散問題。
# 既存N03や共通パッケージには依存しません。
const REQUIRED_METHODS = ("forward_euler", "midpoint_rk2", "classical_rk4", "backward_euler")

"""
    method_ids(include_cn=false)

比較する手法の固定IDをタプルで返す。
`include_cn=true` のときだけ、必須4手法に任意のCrank–Nicolsonを加える。
このIDを計算の選択・図の凡例・summary.tomlで共通に使う。
"""
method_ids(include_cn=false) = include_cn ? (REQUIRED_METHODS..., "crank_nicolson") : REQUIRED_METHODS

"""
    finite_real(x)

`x` が有限な実数なら `true` を返す。数値の指定として `Bool` は受け付けない。
"""
finite_real(x) = x isa Real && !(x isa Bool) && isfinite(x)

"""
    valid_state(y)

状態が有限実数スカラー、または非空で全要素が有限な `Vector{Float64}` か判定する。
"""
valid_state(y) = finite_real(y) || (y isa Vector{Float64} && !isempty(y) && all(isfinite,y))

"""
    state_copy(y)

状態の値を返す。ベクトルなら独立したコピーを作り、入力の配列を変更しない。
"""
state_copy(y::Real) = y
state_copy(y::Vector{Float64}) = copy(y)

"""
    validate_step(y, t, dt)

1ステップの計算に必要な状態・時刻・時間刻みを検証する。

# 引数

- `y`: 有限な実数スカラー（`Bool` を除く）、または空でない有限な `Vector{Float64}`。
- `t`: 有限な時刻。
- `dt`: 有限な正の時間刻み。

不正な入力は `ArgumentError` で拒否する。入力は変更しない。

# 返り値

条件を満たす場合は `nothing`。
"""
function validate_step(y,t,dt)
    valid_state(y) || throw(ArgumentError("状態は有限実数または非空の有限Vector{Float64}が必要です"))
    finite_real(t) || throw(ArgumentError("時刻は有限実数が必要です"))
    finite_real(dt) && dt>0 || throw(ArgumentError("dtは有限正の実数が必要です"))
    nothing
end

"""
    rhs(f, y, t)

右辺 `f(y, t)` を評価し、状態と対応する形の有限な値を返す。
配列の場合は `f` へコピーを渡し、返り値もコピーする。
これにより、右辺が受け取った配列を書き換えても呼出し元の状態を保つ。

`y` と `t` は呼出し元で検証済みとする。
右辺の型・長さが合わなければ `ArgumentError`、非有限値ならエラーで停止する。
"""
function rhs(f,y,t)
    # 評価後は、スカラー／ベクトルの種類とベクトル長が合うか確認する。
    value = f(state_copy(y),t)

    if y isa Vector
        value isa Vector{Float64} && length(value)==length(y) ||
            throw(ArgumentError("RHSは状態と同じ長さのVector{Float64}を返してください"))
    else
        value isa Real || throw(ArgumentError("スカラー状態のRHSは実数を返してください"))
    end

    # 正しい形でもNaNやInfは計算失敗とし、後続の更新へ流さない。
    valid_state(value) || error("RHSが非有限値を返しました")
    state_copy(value)
end

"""
    forward_euler_step(f, y, t, dt)

前進Euler法で `y` を1ステップ進めた新しい状態を返す。
現在の傾きだけを使う1次の陽解法で、RHS評価は1回。

引数の条件は `rk2_step` と同じ。入力状態を変更せず、配列は独立して返す。
"""
function forward_euler_step(f,y,t,dt)
    validate_step(y,t,dt)

    # 現在の時刻における傾きにdtを掛け、旧状態へ加える。
    y + dt*rhs(f,y,t)
end

"""
    classical_rk4_step(f, y, t, dt)

古典的RK4で `y` を1ステップ進めた新しい状態を返す。
区間の始点・中点・終点の傾きを組み合わせる4次の陽解法で、RHS評価は4回。

引数の条件は `rk2_step` と同じ。入力状態を変更せず、配列は独立して返す。
"""
function classical_rk4_step(f,y,t,dt)
    validate_step(y,t,dt)

    # k1は始点、k2とk3は中点、k4は終点に対応する傾き。
    k1 = rhs(f,y,t)
    k2 = rhs(f,y+dt/2*k1,t+dt/2)
    k3 = rhs(f,y+dt/2*k2,t+dt/2)
    k4 = rhs(f,y+dt*k3,t+dt)

    # 傾きを1:2:2:1の重みで平均し、1ステップ分の変化量にする。
    y + dt/6*(k1+2k2+2k3+k4)
end

"""
    validate_linear(A, y)

線形問題 `y′ = A*y` の係数が状態に対応することを確認する。
スカラー状態なら有限実数、ベクトル状態なら同じ次元の有限実数正方行列を求める。
`y` は検証済みとし、不適合な `A` は `ArgumentError` で拒否する。
条件を満たす場合は `nothing` を返し、入力は変更しない。
"""
function validate_linear(A,y)
    valid = if y isa Real
        finite_real(A)
    else
        A isa AbstractMatrix{<:Real} && size(A)==(length(y),length(y)) && all(isfinite,A)
    end
    valid || throw(ArgumentError("線形係数Aは状態に対応する有限実数または正方行列が必要です"))
    nothing
end

"""
    prepare_linear_step(A, y, t, dt; theta=1.0, factorizer=lu, linear_solver=...)

線形問題の陰解法で、同じ時間刻みを繰り返し使う準備をする。

# 引数

- `A`: `y′ = A*y` の有限実数係数、または状態と同じ次元の正方行列。
- `y`, `t`, `dt`: 検証に使う初期状態・時刻・正の時間刻み。
- `theta`: 後退Eulerは `1.0`、Crank–Nicolsonは `0.5`。
- `factorizer`, `linear_solver`: 行列分解と求解を行う関数。通常は既定値を使う。

# 返り値

旧状態と時刻を受け取る関数 `step(old, time)`。
この関数は入力を変更せず、次時刻の状態を返す。
準備時に分解した `I - theta*dt*A` を保持し、呼出しのたびに再利用する。
逆行列は作らない。
"""
function prepare_linear_step(A,y,t,dt; theta=1.0, factorizer=lu, linear_solver=(F,b)->F\b)
    validate_step(y,t,dt)
    validate_linear(A,y)
    theta in (0.5,1.0) || throw(ArgumentError("thetaは0.5または1が必要です"))

    # スカラー問題も1×1行列として扱い、共通の分解・求解を使う。
    matrix = y isa Real ? reshape([A],1,1) : Matrix(A)
    F = factorizer(I-theta*dt*matrix)

    # 内側の関数が分解Fとdtを覚えておくため、次のステップでも使える。
    function step(old,time)
        validate_step(old,time,dt)

        # 後退Eulerの右辺は旧状態。CNでは旧時刻の寄与も加える。
        b = theta==1.0 ? state_copy(old) : old+(1-theta)*dt*(A*old)
        solution = linear_solver(F, old isa Real ? [b] : b)

        # 元の状態と同じ種類に戻す。ベクトルなら独立した配列を返す。
        old isa Real ? only(solution) : copy(solution)
    end
end

"""
    backward_euler_step(A, y, t, dt)

線形問題 `y′ = A*y` の後退Euler法を1ステップだけ実行し、新しい状態を返す。
`(I - dt*A)*y_new = y_old` を解く1次の陰解法。
引数の条件は `prepare_linear_step` と同じで、入力状態は変更しない。

この単発呼出しでは毎回分解する。
複数ステップの `integrate` では準備を1回にまとめ、分解を再利用する。
"""
function backward_euler_step(A,y,t,dt)
    prepare_linear_step(A,y,t,dt)(y,t)
end

"""
    crank_nicolson_step(A, y, t, dt)

線形問題のCrank–Nicolson法を1ステップだけ実行し、新しい状態を返す。
新旧時刻の右辺を半分ずつ用いる2次の陰解法で、任意比較として扱う。
`(I - dt*A/2)*y_new = (I + dt*A/2)*y_old` を解く。
引数の条件は `prepare_linear_step` と同じで、入力状態は変更しない。
後退Eulerと同様に、`integrate` では分解を再利用する。
"""
function crank_nicolson_step(A,y,t,dt)
    prepare_linear_step(A,y,t,dt;theta=0.5)(y,t)
end

"""
    diffusion_problem()

実験Bの格子・拡散行列・初期値・半離散固有値をまとめた名前付きタプルを返す。

領域は `[0, 2]`、拡散係数は0.1、両端を含む41点、格子幅は0.05。
両端Dirichlet 0を既知として除き、内部39点の状態に空間2次中心差分を適用する。
初期値は低周波モード1に、高周波モード39を振幅0.01で加えたもの。
返り値の固有値は、この行列に対応する半離散厳密解を作るために使う。
"""
function diffusion_problem()
    # 全手法で同じ格子・行列を使い、時間積分の違いだけを比較する。
    D = 0.1; dx = 0.05
    x = collect(0:40).*dx

    # 内部点の2階差分は、下側・対角・上側の3本の対角成分で表せる。
    A = Tridiagonal(fill(D/dx^2,38),fill(-2D/dx^2,39),fill(D/dx^2,38))

    # 内部点j=1,…,39で、低周波の分布に小さい高周波成分を重ねる。
    phi1 = sin.(pi*(1:39)/40)
    phi39 = sin.(39pi*(1:39)/40)

    # 連続PDEの固有値ではなく、上の差分行列の固有値を使う。
    lambda(k) = -4D/dx^2*sin(k*pi/80)^2

    (;D,dx,x,A,phi1,phi39,initial=phi1+0.01phi39,lambda1=lambda(1),lambda39=lambda(39))
end
