using Test
include("wide.jl")
include("strict.jl")

@testset "design comparisons" begin
    # 広い入口でも，流出と固定両端の結果は旧入口と一致する。
    u = [2.,3.,4.]
    @test WideBoundary.apply!(u; left_value=5.) === u
    @test u == [5.,3.,3.]
    u = [2.,3.,4.]
    @test WideBoundary.apply!(u; mode=:fixed, left_value=5., right_value=7.) === u
    @test u == [5.,3.,7.]
    # 境界の種類と指定値の組合せに応じた拒否を確かめる。
    for options in ((;mode=:fixed), (;right_value=7.), (;mode=:unknown))
        u = [2.,3.,4.]
        @test_throws ArgumentError WideBoundary.apply!(u; options...)
        @test u == [2.,3.,4.]
    end
    u = [8.,9.]
    # strictの2点拒否と3点受理を確かめる。
    @test_throws ArgumentError StrictBoundary.apply_boundary!(u; left_value=2.)
    @test u == [8.,9.]
    u = [2.,3.,4.]
    @test StrictBoundary.apply_boundary!(u; left_value=5.) === u
    @test u == [5.,3.,3.]
end
