# F02の配列演習。共通の入力検証 → 平均 → 平均からの偏差の順に読みます。
# mean_temperature と temperature_anomaly のTODOを実装します。入力配列は保持します。

module F02JuliaArraysAndTests

export mean_temperature, temperature_anomaly

"""
    validate_temperatures(values::AbstractVector{<:Real})

温度配列が空でなく、すべて有限値であることを確認する。

# 引数

- `values`: 実数ベクトル。検証だけを行い、要素を変更しない。

空配列、NaN、Infを含む配列は `ArgumentError` で拒否する。

# 返り値

条件を満たす場合は `nothing`。
"""
function validate_temperatures(values::AbstractVector{<:Real})
    isempty(values) && throw(ArgumentError("温度の配列を空にはできません"))
    all(isfinite, values) || throw(
        ArgumentError("温度はすべて有限値にしてください。NaNとInfを取り除いてください"),
    )
    nothing
end

"""
    mean_temperature(values::AbstractVector{<:Real})

温度配列の算術平均を求める。

# 引数

- `values`: 空でない有限な実数ベクトル。要素を変更しない。

入力が不正なら共通検証が `ArgumentError` を投げる。

# 返り値

実装後は全要素の算術平均。整数入力でも小数になり得る。

# 受講生のToDo

合計を要素数で割り、平均を返す。配布時は `未実装 F02: mean_temperature` で停止する。
"""
function mean_temperature(values::AbstractVector{<:Real})
    validate_temperatures(values)
    # TODO(F02): すべての値の合計を求め、要素数で割る。
    error("未実装 F02: mean_temperature")
end

"""
    temperature_anomaly(values::AbstractVector{<:Real})

各温度から平均を引いた偏差の配列を作る。

# 引数

- `values`: 空でない有限な実数ベクトル。入力を変更しない。

入力が不正なら共通検証が `ArgumentError` を投げる。

# 返り値

実装後は入力と対応する偏差を持つ新しい配列。整数入力でも偏差は小数になり得る。

# 受講生のToDo

平均を求め、各要素の偏差を新しい配列で返す。配布時は `未実装 F02: temperature_anomaly` で停止する。
"""
function temperature_anomaly(values::AbstractVector{<:Real})
    validate_temperatures(values)
    # TODO(F02): `value - mean_temperature(values)`を要素とする新しい配列を返す。
    error("未実装 F02: temperature_anomaly")
end

end


if abspath(PROGRAM_FILE) == @__FILE__
    sample = [18.0, 20.0, 22.0]
    println("平均 = ", F02JuliaArraysAndTests.mean_temperature(sample))
    println("偏差 = ", F02JuliaArraysAndTests.temperature_anomaly(sample))
end
