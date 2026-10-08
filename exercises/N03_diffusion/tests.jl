# 必須テストは、数値更新・境界条件・熱量・標準計算を順に確認します。
# 末尾の「自作テスト」では、確認したい性質と許容誤差を自分で決めます。

using Test

# 単独実行時は読み込み、Pkg.test などで読込済みなら既存のモジュールを使う。
if !isdefined(Main,:N03Diffusion)
    include(joinpath(@__DIR__, "run.jl"))
end
import .N03Diffusion

@testset "N03 必須テスト" begin
    @testset "小配列1ステップと旧配列不変" begin
        # 非対称な配列を用意し、更新前の値を別の配列に保存する。
        old = [1., 2., 4., 3.]
        saved = copy(old)
        new = similar(old)

        # 新配列の値と、旧配列が変更されていないことを確認する。
        N03Diffusion.diffusion_step!(new, old, 0.25, 1., 1.; boundary=:insulated)
        @test new == [1.5, 2.25, 3.25, 3.5]
        @test old == saved
    end

    @testset "定数場と固定端" begin
        # 各境界条件と両立する定数場は、1ステップ後も変化しない。
        for (boundary, value) in ((:insulated, 0.7), (:fixed, 0.))
            old = fill(value, 5)
            new = similar(old)
            N03Diffusion.diffusion_step!(new, old, 0.2, 1., 1.; boundary)
            @test new ≈ old
        end

        # 既定の固定温度条件では、更新後の両端が0になる。
        old = [1., 2., 4., 3.]
        new = similar(old)
        N03Diffusion.diffusion_step!(new, old, 0.25, 1., 1.)
        @test new == [0., 2.25, 3.25, 0.]
    end

    @testset "台形則の熱量" begin
        # 断熱条件の1ステップ更新前後で、端点の重みを含む熱量を比較する。
        @test N03Diffusion.thermal_content([1., 2., 4., 3.], 1.) == 8.
        @test N03Diffusion.thermal_content([1.5, 2.25, 3.25, 3.5], 1.) == 8.
    end

    @testset "標準計算の安定性" begin
        for boundary in (:fixed, :insulated)
            # 各境界条件で、既定の終了時刻まで計算する。
            r = N03Diffusion.simulate(; boundary)

            # 終了時刻・安定条件・有限性・温度範囲を確認する。
            @test r.steps * r.dt ≈ r.t_final == 1.
            @test 0 < r.fo <= 0.5
            @test all(isfinite, r.u)
            @test -1e-13 <= r.minimum <= r.maximum <= 1 + 1e-13
        end
    end
end

@testset "N03 自作テスト" begin
    # TODO(自作): 非対称な非定数配列を複数ステップ進め、断熱の台形則熱量保存を確認する。
    # 初期熱量を独立に計算し、点数・ステップ数・値の大きさから許容誤差を決める。
    @test false

    # TODO(自作): 固定温度の正弦解析解と41・81・161点の数値解を同時刻・同じfoで比較する。
    # 誤差の減少と観測次数を確認し、指標・許容範囲・検出限界を説明する。
    @test false
end
