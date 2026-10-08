# F01の課題テスト。提供済みの必須確認 → 受講生が追加する自作テストの順に読みます。
# 必須テストの期待値は保持し、自作テストのTODOを編集します。

using Test

if !isdefined(Main, :F01FirstPullRequest)
    include(joinpath(@__DIR__, "run.jl"))
end

# 提供済みの期待値を使い、実装と独立した基準で確認する。
@testset "F01 必須テスト（配布済み）" begin
    @test F01FirstPullRequest.student_greeting("  Thermofluid  ") == "Hello, Thermofluid!"
    @test_throws ArgumentError F01FirstPullRequest.student_greeting("   ")
end

# 必須テストが扱わない条件を選び、期待値を自分で決める。
@testset "F01 自作テスト" begin
    # TODO(自作): 必須テストとは異なる名前を一つ選び、期待する完全な挨拶文字列を自分で書く。
    @test false
end
