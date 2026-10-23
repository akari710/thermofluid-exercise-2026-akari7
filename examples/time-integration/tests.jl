# 時間積分比較のテスト。読む順: 1ステップ → 入力検証 → 実験A/B → 計算回数。
# RK2が未実装なら最初の検査で停止する。仮の値で通さず、run.jlのTODOを実装する。
# 期待値は更新式や固有モードの増幅率から用意し、結果同士の一致だけでは判定しない。

using Test, LinearAlgebra, TOML

# run.jlは読み込みだけでは実験・保存を行わない。テストから関数を呼んで確認する。
include(joinpath(@__DIR__, "run.jl"))
using .TimeIntegrationExample
const E = TimeIntegrationExample

@testset "RK2: midpoint, time argument, independent result" begin
    # 時刻にも依存する右辺を使い、中点の時刻が計算に反映されるか確認する。
    @test E.rk2_step((y,t)->y+t, 1.0, 0.0, 0.2) ≈ 1.24 atol=1e-13 rtol=1e-13

    # ベクトルの既知の更新値を確認し、旧配列が書き換わっていないか調べる。
    y = [1.0, 2.0]
    z = E.rk2_step((y,t)->[-2y[1]+y[2]+t, -y[2]+2t], y, 0.3, 0.2)
    @test z ≈ [1.04, 1.788] atol=1e-13 rtol=1e-13
    @test y == [1.0, 2.0]
    @test E.rk2_step((y,t)->y+t, 1, 0, 0.2) ≈ 1.24

    # 値が等しいだけでは独立性は分からない。返り値を書き換えて確かめる。
    @test z !== y
    z[1] = 100
    @test y == [1.0, 2.0]
end

@testset "provided one-step methods" begin
    # 同じスカラー問題を、各手法の更新式から求めた期待値と比較する。
    @test E.forward_euler_step((y,t)->-y, 2.0, 0.0, 0.1) ≈ 1.8
    @test E.classical_rk4_step((y,t)->-y, 2.0, 0.0, 0.1) ≈ 1.809675
    @test E.backward_euler_step(-1.0, 2.0, 0.0, 0.1) ≈ 2/1.1
    @test E.crank_nicolson_step(-1.0, 2.0, 0.0, 0.1) ≈ 2*0.95/1.05

    # 陰解法は行列を使う問題でも、線形方程式の解と入力非破壊を確認する。
    y = [1.0, 2.0]; A = [-2.0 1.0; 0.0 -1.0]
    for (step, expected) in ((E.backward_euler_step, (I-0.2A)\y),
                             (E.crank_nicolson_step, (I-0.1A)\((I+0.1A)*y)))
        z = step(A,y,0.3,0.2)
        @test z ≈ expected
        @test z !== y
        @test y == [1.0,2.0]
    end
end

@testset "input validation before RHS or state changes" begin
    # 不正入力では、右辺の呼出しにも状態の更新にも進まないことを確かめる。
    calls = Ref(0)
    f(y,t) = (calls[] += 1; -y)
    for step in (E.forward_euler_step, E.rk2_step, E.classical_rk4_step)
        for bad in (NaN, Inf, Float64[], [1.0,NaN], [1,2], "state")
            @test_throws ArgumentError step(f,bad,0.0,0.1)
        end
        for (t,dt) in ((NaN,0.1),(Inf,0.1),(0.0,0.0),(0.0,-0.1),(0.0,NaN),(0.0,Inf))
            y = [1.0,2.0]
            @test_throws ArgumentError step(f,y,t,dt)
            @test y == [1.0,2.0]
        end
    end
    @test calls[] == 0

    # 陰解法では、係数Aの大きさ・形・有限性も検証の対象になる。
    for step in (E.backward_euler_step, E.crank_nicolson_step)
        for A in (zeros(3,3), [NaN 0.0; 0.0 1.0], [1.0 2.0], "matrix")
            y = [1.0,2.0]
            @test_throws ArgumentError step(A,y,0.0,0.1)
            @test y == [1.0,2.0]
        end
        @test_throws ArgumentError step(-1.0,NaN,0.0,0.1)
        @test_throws ArgumentError step(-1.0,1.0,0.0,-0.1)
    end

    # 共通の時間ループの入力も、実際に計算を始める前に拒否する。
    for kwargs in ((;N=0),(;N=true),(;N=2.5),(;T=0.0),(;T=Inf),(;method="unknown"))
        @test_throws ArgumentError E.integrate(f,1.0;kwargs...)
    end
    @test calls[] == 0
end

@testset "A: time convergence and actual work counts" begin
    # 4種類のNから3区間の観測次数を求め、すべての区間が範囲内か調べる。
    data = E.experiments(;include_cn=true)
    ranges = Dict("forward_euler"=>(0.85,1.15), "midpoint_rk2"=>(1.85,2.15),
                  "classical_rk4"=>(3.8,4.2), "backward_euler"=>(0.85,1.15),
                  "crank_nicolson"=>(1.85,2.15))

    # 陽解法は1ステップの段数、陰解法は分解1回とN回の求解が基準。
    stages = Dict("forward_euler"=>1,"midpoint_rk2"=>2,"classical_rk4"=>4)
    for row in data["time_convergence"]
        id = row["method"]; lo,hi = ranges[id]

        # 誤差・次数に加え、全手法で同じdtと終端時刻を使うことも確認する。
        @test length(row["orders"]) == 3
        @test all(p->lo<p<hi,row["orders"])
        @test all(e->isfinite(e)&&e>0,row["errors"])
        @test row["N"] == [8,16,32,64]
        @test row["dt"] == 1.0 ./ row["N"]
        @test row["final_time"] == ones(4)

        for (N,cost) in zip(row["N"],row["cost"])
            explicit = haskey(stages,id)
            @test cost["rhs_evaluations"] == (explicit ? stages[id]*N : 0)
            @test cost["factorizations"] == (explicit ? 0 : 1)
            @test cost["linear_solves"] == (explicit ? 0 : N)
        end
    end
end

@testset "B: semidiscrete eigenmodes and finite amplification" begin
    # 境界を除いた39点の行列が、狙った固有モードを持つことを確認する。
    grid = E.diffusion_problem()
    @test length(grid.x) == 41
    @test size(grid.A) == (39,39)
    @test grid.dx ≈ 0.05

    # 固有モードでは行列問題をz=lambda*dtのスカラー増幅率で調べられる。
    amplification(id,z) = id=="forward_euler" ? 1+z : id=="midpoint_rk2" ? 1+z+z^2/2 :
        id=="classical_rk4" ? 1+z+z^2/2+z^3/6+z^4/24 : id=="backward_euler" ? 1/(1-z) : (1+z/2)/(1-z/2)
    for k in (1,39)
        phi = sin.(k*pi*(1:39)/40)
        lambda = -4*0.1/0.05^2*sin(k*pi/80)^2
        @test grid.A*phi ≈ lambda*phi atol=3e-12

        # 1ステップはg倍、Nステップはg^N倍になることを独立に確認する。
        for id in E.method_ids(true), N in (128,4)
            result = E.integrate((y,t)->grid.A*y,phi;A=grid.A,T=0.1,N,method=id)
            g = amplification(id,lambda*0.1/N)
            @test result.states[2] ≈ g*phi atol=3e-13
            @test result.states[end] ≈ g^N*phi atol=5e-11 rtol=2e-11
            @test phi ≈ sin.(k*pi*(1:39)/40)
        end
    end

    # 2つのモードを重ねた実験結果も、条件・参照解・係数履歴まで確認する。
    data = E.experiments(;include_cn=true)
    for row in data["diffusion"]
        N = row["N"]; id = row["method"]; dt = row["dt"]
        @test N in (128,4)
        @test dt == 0.1/N
        @test row["Fo"] ≈ (N==128 ? 0.03125 : 1.0)
        @test row["final_time"] == 0.1
        @test row["times"] == collect(0:N).*dt
        @test row["final_state"][[1,end]] == [0.0,0.0]
        @test row["initial_state"] ≈ vcat(0.0,grid.initial,0.0)

        # 連続PDEの解ではなく、同じ差分行列の固有値から参照値を作る。
        reference = vcat(0.0,exp(grid.lambda1*0.1)*grid.phi1+0.01exp(grid.lambda39*0.1)*grid.phi39,0.0)
        @test row["reference_state"] ≈ reference
        @test row["max_abs_error"] ≈ maximum(abs.(row["final_state"]-reference))

        # モード39の全時刻の係数を、初期振幅0.01と増幅率から予測する。
        g = amplification(id,grid.lambda39*dt)
        @test row["high_mode_coefficients"] ≈ 0.01 .* g.^(0:N) atol=2e-11 rtol=2e-10
        @test all(isfinite,row["final_state"])

        # 粗いdtでの有限な増幅は意図した結果。CNは減衰しつつ符号が交代する。
        if N==4
            if id in ("forward_euler","midpoint_rk2","classical_rk4")
                @test abs(last(row["high_mode_coefficients"])) > 0.01
            elseif id=="crank_nicolson"
                @test all(row["high_mode_coefficients"][1:end-1].*row["high_mode_coefficients"][2:end] .< 0)
                @test abs(last(row["high_mode_coefficients"])) < 0.01
            end
        end
    end
end

@testset "work counters measure calls, not diagnostics" begin
    # 検査側でも実呼出しを数え、記録された回数が作られた値でないか確かめる。
    for id in E.method_ids(true)
        evaluations = Ref(0); factorizations = Ref(0); solves = Ref(0)
        f(y,t) = (evaluations[]+=1; -y)
        factorizer = M->(factorizations[]+=1; lu(M))
        solver = (F,b)->(solves[]+=1; F\b)

        r = E.integrate(f,1.0;method=id,N=8,T=1.0,A=-1.0,factorizer,linear_solver=solver)
        @test r.cost == Dict("rhs_evaluations"=>evaluations[],"factorizations"=>factorizations[],"linear_solves"=>solves[])
    end
end

@testset "independent input even when RHS writes its argument" begin
    # 右辺が受け取った配列を書き換えても、元の初期状態には影響させない。
    y = [1.0,2.0]
    for step in (E.forward_euler_step,E.rk2_step,E.classical_rk4_step)
        z = step((u,t)->(u[1]=3.0; u),y,0.0,0.1)
        @test y == [1.0,2.0]
        @test z !== y
    end

    # スカラーの初期値が整数でも、浮動小数の計算結果を履歴に保存できる。
    @test E.integrate((y,t)->-y,1;method="midpoint_rk2").states[end] ≈ E.integrate((y,t)->-y,1.0;method="midpoint_rk2").states[end]
end
