# 必須テストは提供済みの期待値で確認します。末尾の自作テストだけを受講生が編集します。
# 入力の準備 → 独立した期待値との比較 → 入力保持・不正入力の確認の順に読みます。

using Test
if !isdefined(Main, :N04AdvectionDiffusion)
    include(joinpath(@__DIR__, "run.jl"))
end
import .N04AdvectionDiffusion
const N04 = N04AdvectionDiffusion
N04.validate_model(N04.SELECTED_MODEL)
@testset "N04 必須テスト" begin
    model = N04.SELECTED_MODEL
    @testset "非対称1ステップ・周期端・旧配列不変" begin
        old = [1., 2., 4., 3.]
        saved = copy(old)
        new = similar(old)
        expected = model == :linear ? [1.26, 1.92, 3.74, 3.08] : [1.46, 1.87, 3.34, 3.33]
        @test N04.advection_diffusion_step!(new, old, 0.1, 1., 0.2) === new
        @test new ≈ expected atol = 1e-14
        @test old == saved
    end
    @testset "定数場と周期積分" begin
        old = fill(1.7, 5)
        new = similar(old)
        N04.advection_diffusion_step!(new, old, 0.1, 1., 0.2)
        @test new ≈ old
        @test N04.conserved_integral([1., 2., 4., 3.], 0.5) == 5.
    end
    @testset "移流・拡散をゼロにした極限" begin
        old = [1., 2., 4., 3.]
        new = similar(old)
        N04.advection_diffusion_step!(new, old, 0.1, 1., 0.2; advection=false)
        @test new ≈ [1.06, 2.02, 3.94, 2.98]
        N04.advection_diffusion_step!(new, old, 0.1, 1., 0.)
        @test new ≈ (model == :linear ? [1.2, 1.9, 3.8, 3.1] : [1.4, 1.85, 3.4, 3.35])
    end
    @testset "合成刻み・最終時刻・標準安定性" begin
        @test N04.stable_timestep(2., 0.5, 0.1) ≈ 1 / 6
        r = N04.simulate()
        @test r.steps * r.dt ≈ r.t_final == 1.
        @test r.max_stability_number <= 0.8 + 1e-12
        @test all(isfinite, r.u) && 1 - 1e-12 <= r.minimum <= r.maximum <= 2 + 1e-12
    end
end
@testset "N04 自作テスト" begin
    # TODO(自作): 非対称な非定数配列を複数ステップ進め、周期積分の保存を確認する。
    # 初期積分を独立に計算し、点数・ステップ数・値の大きさから許容誤差を説明する。
    # TODO(N04): 複数ステップ保存
    @testset "複数ステップ保存" begin
        # 非対称・非定数・非負の初期配列
        old = [1.0, 1.4, 2.0, 1.2, 1.7]
        new = similar(old)

        dx = 0.4
        diffusivity = 0.1
        steps = 100

        # 初期最大値を速度上限として，安定な時間刻みを設定する。
        # 線形では標準速度1.0，非線形では初期配列の最大値を使う。
        max_speed = N04.SELECTED_MODEL == :linear ? 1.0 : maximum(old)
        dt = N04.stable_timestep(max_speed, dx, diffusivity)

        # conserved_integralを使わず，定義から初期積分を計算する。
        initial_integral = dx * sum(old)

        # 5点，100ステップ，値が1～2程度であることを踏まえる
        # 丸め誤差には余裕を持たせた絶対許容誤差とする。
        tolerance = 1e-12

        for step in 1:steps
            N04.advection_diffusion_step!(
                new, old, dt, dx, diffusivity
            )


            # 各ステップで初期積分との一致を確認する。
            @test isapprox(
                dx * sum(new),
                initial_integral;
                atol=tolerance,
                rtol=0.0,
            )

            old, new = new, old
        end
    end
    # TODO(自作): 選択モデルの解析解と40・80・160点の誤差・観測次数を確認する。
    # 誤差の指標、同じ最終時刻、0.8〜1.2の次数範囲、検出できない誤りも説明する。
    # TODO(N04): 格子収束
    @testset "解析解と格子収束" begin
        errors = Float64[]

        for nx in (40, 80, 160)
            r = N04.simulate(;
                nx=nx,
                diffusivity=0.1,
                safety=0.8,
                t_final=1.0,
                initial=:mode,
            )

            # 全格子で同じ最終時刻に到達したことを確認する。
            @test r.steps * r.dt ≈ 1.0

            # 各格子点の解析解を計算する。
            exact = [
                N04.analytic_solution(
                    x, r.t_final;
                    model=N04.SELECTED_MODEL,
                    diffusivity=0.1,
                    speed=1.0,
                )
                for x in r.x
            ]

            # 誤差指標は最大絶対誤差。
            error = maximum(abs.(r.u .- exact))
            push!(errors, error)
        end

        # 格子を細かくすると誤差が減少することを確認する。
        @test all(isfinite, errors)
        @test 0 < errors[3] < errors[2] < errors[1]

        # 40→80点，80→160点の観測次数
        p1 = log2(errors[1] / errors[2])
        p2 = log2(errors[2] / errors[3])

        @test 0.8 < p1 < 1.2
        @test 0.8 < p2 < 1.2
    end
end
