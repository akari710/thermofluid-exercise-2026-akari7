include("common.jl")
"""
流入・流出の旧入口。
left_valueを共通側へ渡し，その戻り値を返す。
"""
function apply_boundary!(u::Vector{Float64}; left_value::Float64=1.0)
    # 省略時も明示時も，同じキーワードで値を渡す。
    return BoundaryCore.inlet_outlet!(u; left_value=left_value)
end
"""
固定両端の旧入口。
左右の指定値を共通側へ渡し，その戻り値を返す。
"""
function fixed_boundary!(u::Vector{Float64};
    left_value::Float64, right_value::Float64)
    return BoundaryCore.fixed_ends!(u;
        left_value=left_value, right_value=right_value)
end
