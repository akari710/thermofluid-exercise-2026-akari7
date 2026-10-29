# 必須テストは提供済みの期待値で確認します。末尾の自作テストだけを受講生が編集します。
# 入力の準備 → 独立した期待値との比較 → 入力保持・不正入力の確認の順に読みます。

using Test
if !isdefined(Main, :N02NonlinearAdvection)
    include(joinpath(@__DIR__, "run.jl"))
end
using .N02NonlinearAdvection
@testset "N02 必須テスト" begin
    @testset "流束と周期添字" begin
        @test burgers_flux(1.) == 0.5
        @test burgers_flux(2.) == 2.
        @test periodic_left_index(1, 4) == 4
        @test periodic_left_index(3, 4) == 2
    end
    @testset "保存形1ステップと旧配列不変" begin
        old = [1., 2., 1., 1.]
        saved = copy(old)
        new = similar(old)
        nonlinear_upwind_step!(new, old, 0.25, 1.; boundary=:periodic)
        @test new == [1., 1.625, 1.375, 1.]
        @test old == saved
    end
    @testset "定数場と境界条件" begin
        for (boundary, value) in ((:fixed, 1.), (:periodic, 0.7))
            old = fill(value, 5)
            new = similar(old)
            nonlinear_upwind_step!(new, old, 0.2, 1.; boundary)
            apply_boundary!(new; boundary)
            @test new ≈ old
        end
        @test apply_boundary!([4., 2., 3., 8.]) == [1., 2., 3., 3.]
    end
    @testset "CFL・有界性・非線形な波形変形" begin
        for (boundary, nx) in ((:fixed, 81), (:periodic, 80))
            r = simulate(; boundary, nx)
            @test r.max_cfl <= 0.5 + 1e-14
            @test 1 - 1e-13 <= minimum(r.u) <= maximum(r.u) <= 2 + 1e-13
            # c=1で同じ最終時刻まで移した矩形との離散L1差。周期では左右をつなぐ。
            shifted = [
                0.5 <=
                    (boundary == :periodic ? mod(x - r.t_final, 2.) : x - r.t_final) <=
                    1.0 ? 2. : 1. for x in r.x
            ]
            @test r.dx * sum(abs.(r.u - shifted)) > 0.2
            @test maximum(r.u[findall(>(1.5), r.x)]) > 1.2
        end
    end
end
@testset "N02 自作テスト" begin
    # TODO(自作): 周期の非定数・非負配列を複数ステップ進め、離散総和保存を確認する。入力・許容誤差・期待値を自分で決める。
    @testset "周期境界での総和保存" begin
        # 非定数で，すべての値が0以上の配列を用意
        u0 = [0.2, 1.0, 0.4, 0.0, 1.8, 0.7, 0.1, 0.5]

        # 時間刻みと空間刻み
        dt = 0.02
        dx = 0.1

        old = copy(u0)
        new = similar(old)

        # 複数ステップ計算する
        for n in 1:20
            nonlinear_upwind_step!(new, old, dt, dx; boundary=:periodic)
            old, new = new, old
        end

        # 計算前後で総和が保存されることを確認
        @test isapprox(
            sum(old), sum(u0);
            atol=1e-10,
            rtol=1e-10
        )
    end

    # TODO(自作): 周期1ステップSについてS(circshift(u,k)) ≈ circshift(S(u),k)を確認する。非対称な非定数配列、k≠0、同じdtとdxを使う。
    @testset "周期境界での平行移動対称性" begin
        # 非対称な非定数配列
        u = [0.2, 1.3, 0.0, 0.7, 2.1, 0.4, 0.9, 0.1]

        # ずらす要素数
        k = 3

        # 両方の計算で同じ値を使う
        dt = 0.02
        dx = 0.1

        # 先にずらしてから1ステップ計算する
        shifted_u = circshift(u, k)
        shifted_then_step = similar(u)

        nonlinear_upwind_step!(
            shifted_then_step, shifted_u, dt, dx;
            boundary=:periodic
        )

        # 先に1ステップ計算してからずらす
        stepped_u = similar(u)

        nonlinear_upwind_step!(
            stepped_u, u, dt, dx;
            boundary=:periodic
        )

        step_then_shifted = circshift(stepped_u, k)

        # 二つの結果が一致することを確認
        @test isapprox(
            shifted_then_step, step_then_shifted;
            atol=1e-10,
            rtol=1e-10
        )
    end

end
