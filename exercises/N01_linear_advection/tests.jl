# N01のテスト枠。必須テストと自作テストの両方を受講生が実装します。
# 各TODOで入力と独立した期待値を選び、配列全体・境界・安定性を確認します。

using Test

if !isdefined(Main, :N01LinearAdvection)
    include(joinpath(@__DIR__, "run.jl"))
end

@testset "N01 必須テスト" begin
    # TODO(必須): 座標、base、plateau、区間端を自分で選び、rectangular_initial_conditionが返す配列全体を手計算した期待値と比較する。
    @test false

    # TODO(必須): 小さなu_old、移流速度、時間刻み、格子幅を選び、upwind_step!後のu_new配列全体を手計算した期待値と比較する。
    @test false

    # TODO(必須): 同じ入力でcentered_step!後のu_new配列全体を手計算した期待値と比較する。風上差分との違いが現れる入力を選ぶ。
    @test false

    # TODO(必須): apply_boundary!が左端と右端を課題の境界条件どおりに書き換えることを、配列全体の期待値で確かめる。
    @test false

    # TODO(必須): simulateで安定な風上差分と不安定な中心差分を同じ条件で計算し、最終値の範囲または振幅から両者の違いを確かめる。条件と判定基準は自分で書く。
    @test false
end

@testset "N01 自作テスト" begin
    # TODO(自作): 必須とは異なる初期条件、CFL、格子数、またはAPIの性質を一つ選び、期待する結果を自分で書く。PNG、TOML、ファイルサイズは対象にしない。
    @test false
end
