module BoundaryCore
# 検証も同じモジュールへ移し，共通側から呼べるようにする。
"""
配列長2点以上と境界値の有限性を確認する。
成功時はnothingを返す。
"""
function validate_boundary_inputs(u, values...)
    length(u) >= 2 || throw(ArgumentError("2点以上が必要です"))
    all(isfinite, values) || throw(ArgumentError("境界値は有限にします"))
    return nothing
end

"""
    inlet_outlet!(u; left_value=1.0)

左端へleft_value，右端へ左隣の値を書き，同じ配列uを返す。
"""
function inlet_outlet!(u::Vector{Float64}; left_value::Float64=1.0)
    validate_boundary_inputs(u, left_value)
    # 抽出前と同じ順序で，左端から右端へ更新する。
    u[1] = left_value
    u[end] = u[end - 1]
    return u
end

"""
    fixed_ends!(u; left_value, right_value)

左右の端へ指定値を書き，内部点を保った同じ配列uを返す。
"""
function fixed_ends!(u::Vector{Float64};
    left_value::Float64, right_value::Float64)
    validate_boundary_inputs(u, left_value, right_value)
    u[1] = left_value
    u[end] = right_value
    return u
end
end
