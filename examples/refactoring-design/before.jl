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
    apply_boundary!(u; left_value=1.0)

左端へleft_value，右端へ左隣の値を書き，同じ配列uを返す。
"""
function apply_boundary!(u::Vector{Float64}; left_value::Float64=1.0)
    # 配列長と境界値を，書込み前に確認する。
    validate_boundary_inputs(u, left_value)
    u[1] = left_value
    # 2点配列では，直前に更新した左端を読む。
    u[end] = u[end - 1]
    return u
end

"""
左右の端へ指定値を書き，内部点を保つ。
戻り値は更新した同じ配列u。
"""
function fixed_boundary!(u::Vector{Float64};
    left_value::Float64, right_value::Float64)
    validate_boundary_inputs(u, left_value, right_value)
    u[1] = left_value
    u[end] = right_value
    return u
end
