using Test
# 実行引数で比較対象を選び，全対象に同じ検査を適用する。
stage = only(ARGS)
include(stage * ".jl")
# 共通側の直接呼出しと旧入口を，同じf・gとして扱う。
f = stage == "common" ? BoundaryCore.inlet_outlet! : apply_boundary!
g = stage == "common" ? BoundaryCore.fixed_ends! : fixed_boundary!

@testset "boundary contract: $stage" begin
    # 既定値・明示値・2点での更新順序を確かめる。
    for (input, options, expected) in (
        ([2.,3.,4.], (;), [1.,3.,3.]),
        ([2.,3.,4.], (;left_value=5.), [5.,3.,3.]),
        ([8.,9.], (;left_value=2.), [2.,2.]),
    )
        u = copy(input)
        # 配列の同一性と，更新後の値をそれぞれ確かめる。
        @test f(u; options...) === u
        @test u == expected
    end
    # 拒否する入力では，例外と元の配列の保持を確かめる。
    for input in (Float64[], [8.], [8.,9.], [2.,3.,4.])
        for value in (1., Inf, NaN)
            if length(input) < 2 || !isfinite(value)
                u = copy(input)
                @test_throws ArgumentError f(u; left_value=value)
                @test isequal(u, input)
            end
        end
        for (left, right) in ((5.,7.), (Inf,7.), (NaN,7.), (5.,Inf), (5.,NaN))
            if length(input) < 2 || !all(isfinite, (left,right))
                u = copy(input)
                @test_throws ArgumentError g(u; left_value=left, right_value=right)
                @test isequal(u, input)
            end
        end
    end
    # 固定両端では，左右それぞれの指定値を確かめる。
    u = [2.,3.,4.]
    @test g(u; left_value=5., right_value=7.) === u
    @test u == [5.,3.,7.]
end
