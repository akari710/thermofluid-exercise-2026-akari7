# F03-F04の課題テスト。提供済みの必須確認 → 受講生が追加する自作テストの順に読みます。
# 必須テストの期待値は保持し、自作テストのTODOを編集します。

using Test

if !isdefined(Main, :F04NumericalDifferentiation)
    include(joinpath(@__DIR__, "run.jl"))
end
using .F04NumericalDifferentiation

# 提供済みの期待値を使い、実装と独立した基準で確認する。
@testset "F03-F04 必須テスト（配布済み）" begin
    @testset "差分値（forward_difference・backward_difference・centered_difference）" begin

        """
            quadratic(x)

        必須テストで使う2次関数を評価する。

        # 引数

        - `x`: テストの評価点。

        # 返り値

        指定点の参照値。
        """
        quadratic(x) = x ^ 2
        @test forward_difference(quadratic, 2.0, 0.5) ≈ 4.5
        @test backward_difference(quadratic, 2.0, 0.5) ≈ 3.5
        @test centered_difference(quadratic, 2.0, 0.5) ≈ 4.0
    end

    @testset "収束次数（convergence_study）" begin

        """
            cubic(x)

        必須テストで使う3次関数を評価する。

        # 引数

        - `x`: テストの評価点。

        # 返り値

        指定点の参照値。
        """
        cubic(x) = x ^ 3

        """
            cubic_derivative(x)

        必須テストで使う3次関数の解析的な微分を評価する。

        # 引数

        - `x`: テストの評価点。

        # 返り値

        指定点の参照値。
        """
        cubic_derivative(x) = 3x ^ 2
        study = convergence_study(cubic, cubic_derivative, 1.0, [0.2, 0.1, 0.05])
        @test all(ratio -> 1.8 <= ratio <= 2.2, study.forward_ratios)
        @test all(ratio -> 1.8 <= ratio <= 2.2, study.backward_ratios)
        @test all(ratio -> 3.9 <= ratio <= 4.1, study.centered_ratios)
    end

    @testset "ベクトル公式の格子収束（verify_vector_identities）" begin
        coarse = verify_vector_identities(9)
        fine = verify_vector_identities(17)
        @test keys(coarse) ==
              (:curl_gradient, :divergence_curl, :product_divergence, :curl_curl)
        for key in keys(coarse)
            @test 0 < fine[key] < coarse[key]
            @test 3.0 <= coarse[key] / fine[key] <= 4.8
        end
    end
end

# 必須テストが扱わない条件を選び、期待値を自分で決める。
@testset "F03-F04 自作テスト" begin
    # TODO(自作): 別の関数、評価点、または入力条件を選び、どの実装ミスを検出するか説明できるテストを一つ書く。F03単独の別テストは作らない。
    @test false
end
