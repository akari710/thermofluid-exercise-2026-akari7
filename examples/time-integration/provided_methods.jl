# 配布実装。既存N03や共通パッケージには依存しない。
const REQUIRED_METHODS = ("forward_euler", "midpoint_rk2", "classical_rk4", "backward_euler")
method_ids(include_cn=false) = include_cn ? (REQUIRED_METHODS..., "crank_nicolson") : REQUIRED_METHODS

finite_real(x) = x isa Real && !(x isa Bool) && isfinite(x)
valid_state(y) = finite_real(y) || (y isa Vector{Float64} && !isempty(y) && all(isfinite,y))
state_copy(y::Real) = y
state_copy(y::Vector{Float64}) = copy(y)

function validate_step(y,t,dt)
    valid_state(y) || throw(ArgumentError("状態は有限実数または非空の有限Vector{Float64}が必要です"))
    finite_real(t) || throw(ArgumentError("時刻は有限実数が必要です"))
    finite_real(dt) && dt>0 || throw(ArgumentError("dtは有限正の実数が必要です"))
    nothing
end

# RHSにも独立な状態を渡す。学生の入力配列を保持する。
function rhs(f,y,t)
    value = f(state_copy(y),t)
    if y isa Vector
        value isa Vector{Float64} && length(value)==length(y) ||
            throw(ArgumentError("RHSは状態と同じ長さのVector{Float64}を返してください"))
    else
        value isa Real || throw(ArgumentError("スカラー状態のRHSは実数を返してください"))
    end
    valid_state(value) || error("RHSが非有限値を返しました")
    state_copy(value)
end

function forward_euler_step(f,y,t,dt)
    validate_step(y,t,dt)
    y + dt*rhs(f,y,t)
end

function classical_rk4_step(f,y,t,dt)
    validate_step(y,t,dt)
    k1 = rhs(f,y,t)
    k2 = rhs(f,y+dt/2*k1,t+dt/2)
    k3 = rhs(f,y+dt/2*k2,t+dt/2)
    k4 = rhs(f,y+dt*k3,t+dt)
    y + dt/6*(k1+2k2+2k3+k4)
end

function validate_linear(A,y)
    valid = if y isa Real
        finite_real(A)
    else
        A isa AbstractMatrix{<:Real} && size(A)==(length(y),length(y)) && all(isfinite,A)
    end
    valid || throw(ArgumentError("線形係数Aは状態に対応する有限実数または正方行列が必要です"))
    nothing
end

# 固定dtの実験ではこの準備を1回だけ呼び、分解を全ステップで再利用する。
function prepare_linear_step(A,y,t,dt; theta=1.0, factorizer=lu, linear_solver=(F,b)->F\b)
    validate_step(y,t,dt)
    validate_linear(A,y)
    theta in (0.5,1.0) || throw(ArgumentError("thetaは0.5または1が必要です"))
    matrix = y isa Real ? reshape([A],1,1) : Matrix(A)
    F = factorizer(I-theta*dt*matrix)
    function step(old,time)
        validate_step(old,time,dt)
        b = theta==1.0 ? state_copy(old) : old+(1-theta)*dt*(A*old)
        solution = linear_solver(F, old isa Real ? [b] : b)
        old isa Real ? only(solution) : copy(solution)
    end
end

function backward_euler_step(A,y,t,dt)
    prepare_linear_step(A,y,t,dt)(y,t)
end

function crank_nicolson_step(A,y,t,dt)
    prepare_linear_step(A,y,t,dt;theta=0.5)(y,t)
end

function diffusion_problem()
    D = 0.1; dx = 0.05
    x = collect(0:40).*dx
    A = Tridiagonal(fill(D/dx^2,38),fill(-2D/dx^2,39),fill(D/dx^2,38))
    phi1 = sin.(pi*(1:39)/40)
    phi39 = sin.(39pi*(1:39)/40)
    lambda(k) = -4D/dx^2*sin(k*pi/80)^2
    (;D,dx,x,A,phi1,phi39,initial=phi1+0.01phi39,lambda1=lambda(1),lambda39=lambda(39))
end
