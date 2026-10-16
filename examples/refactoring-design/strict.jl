module StrictBoundary
include("common.jl")
"""
3点以上へ受理条件を変えた比較例。
共通側の更新結果を返す。
"""
function apply_boundary!(u::Vector{Float64}; left_value::Float64=1.0)
    # 従来有効だった2点入力が，共通側を呼ぶ前に拒否される。
    length(u) >= 3 || throw(ArgumentError("3点以上が必要です"))
    return BoundaryCore.inlet_outlet!(u; left_value=left_value)
end
end
