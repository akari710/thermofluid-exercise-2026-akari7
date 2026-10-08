# F02の課題テスト。提供済みの必須確認 → 受講生が追加する自作テストの順に読みます。
# 必須テストの期待値は保持し、自作テストのTODOを編集します。

using Test

if !isdefined(Main, :F02JuliaArraysAndTests)
    include(joinpath(@__DIR__, "run.jl"))
end

# 提供済みの期待値を使い、実装と独立した基準で確認する。
@testset "F02 必須テスト（配布済み）" begin
    values = [5.0, 7.0, 12.0]
    original = copy(values)
    anomalies = F02JuliaArraysAndTests.temperature_anomaly(values)

    # この具体例の平均と偏差は厳密に表せる値なので、==で比較する。
    @test F02JuliaArraysAndTests.mean_temperature(values) == 8.0
    @test anomalies == [-3.0, -1.0, 4.0]

    # originalは呼び出し前にcopyした値。単なる代入では変更を見逃す。
    @test values == original
    @test_throws ArgumentError F02JuliaArraysAndTests.mean_temperature(Float64[])
end

# 必須テストが扱わない条件を選び、期待値を自分で決める。
@testset "F02 自作テスト" begin
    # TODO(自作): 戻り値の型、別の数学的性質、または必須とは異なる不正入力から一つ選び、入力と期待値を自分で書く。
    # 整数入力でも平均や偏差は小数になり得る。入力と同じ型かではなく、計算結果を格納できるかを考える。
    # Float32を使う場合は、値の一致とtypeof／eltypeによる型の確認を区別する。
    @test false
end
