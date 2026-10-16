module WideBoundary
include("common.jl")
"""
    apply!(u; mode=:outlet, left_value=1.0, right_value=nothing)

modeで境界処理を選び，更新した同じ配列uを返す。
:outletでは左値を使い，:fixedでは左右の値を使う。
"""
function apply!(u::Vector{Float64}; mode=:outlet,
    left_value::Float64=1.0, right_value=nothing)
    # 境界の選択と，有効な指定の組合せをここで検証する。
    if mode == :outlet
        isnothing(right_value) || throw(ArgumentError("流出で指定する境界値は左値だけです"))
        return BoundaryCore.inlet_outlet!(u; left_value=left_value)
    elseif mode == :fixed
        right_value isa Float64 || throw(ArgumentError("固定では右値が必要です"))
        return BoundaryCore.fixed_ends!(u;
            left_value=left_value, right_value=right_value)
    else
        throw(ArgumentError("未知の境界です"))
    end
end
end
