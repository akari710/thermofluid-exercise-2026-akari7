# F03の解析的な場と微分結果を提供します。場の定義 → 点の検証 → 各微分量の順に読みます。
# F04ではこれらを独立した参照値として使います。このファイルの数値式は編集しません。

module F03VectorCalculus

export scalar_field,
    vector_field,
    gradient_scalar,
    curl_vector,
    laplacian_scalar,
    divergence_vector,
    gradient_divergence_vector,
    laplacian_vector

"""
    scalar_field(point)

教材で指定した3次元のスカラー場を評価する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

# 返り値

指定点におけるスカラー場の値。座標の検証はこの関数では行わない。
"""
scalar_field(point) = sin(point[1]) * cos(point[2]) * exp(point[3])

"""
    vector_field(point)

教材で指定した3次元のベクトル場を評価する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

# 返り値

指定点の3成分を持つタプル。座標の検証はこの関数では行わない。
"""
vector_field(point) = (
    sin(point[2]) * exp(2point[3]) + exp(point[1]),
    sin(point[3]) * exp(2point[1]) + exp(point[2]),
    sin(point[1]) * exp(2point[2]) + exp(point[3]),
)

"""
    validate_point(point)

解析的な微分量を評価する座標を検証する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

3要素のタプルでない場合、または有限な実数でない座標は `ArgumentError` で拒否する。

# 返り値

条件を満たす場合は `nothing`。
"""
function validate_point(point)
    point isa Tuple && length(point) == 3 ||
        throw(ArgumentError("点は3要素のタプルで指定してください"))
    all(value -> value isa Real && isfinite(value), point) ||
        throw(ArgumentError("点の座標は有限な実数にしてください"))
    nothing
end

"""
    gradient_scalar(point)

スカラー場の解析的な勾配を評価する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

評価前に `validate_point` で有限な実数座標か確認する。不正な入力は `ArgumentError`。

# 返り値

x・y・z方向の偏微分を持つ3要素のタプル。
"""
function gradient_scalar(point)
    validate_point(point)
    (
        cos(point[1]) * cos(point[2]) * exp(point[3]),
        -sin(point[1]) * sin(point[2]) * exp(point[3]),
        sin(point[1]) * cos(point[2]) * exp(point[3]),
    )
end

"""
    curl_vector(point)

ベクトル場の解析的な回転を評価する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

評価前に `validate_point` で有限な実数座標か確認する。不正な入力は `ArgumentError`。

# 返り値

x・y・z方向の回転成分を持つ3要素のタプル。
"""
function curl_vector(point)
    validate_point(point)
    (
        2sin(point[1]) * exp(2point[2]) - cos(point[3]) * exp(2point[1]),
        2sin(point[2]) * exp(2point[3]) - cos(point[1]) * exp(2point[2]),
        2sin(point[3]) * exp(2point[1]) - cos(point[2]) * exp(2point[3]),
    )
end

"""
    laplacian_scalar(point)

スカラー場の解析的なラプラシアンを評価する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

評価前に `validate_point` で有限な実数座標か確認する。不正な入力は `ArgumentError`。

# 返り値

指定点のラプラシアンのスカラー値。
"""
function laplacian_scalar(point)
    validate_point(point)
    -scalar_field(point)
end

"""
    divergence_vector(point)

ベクトル場の解析的な発散を評価する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

評価前に `validate_point` で有限な実数座標か確認する。不正な入力は `ArgumentError`。

# 返り値

指定点の発散のスカラー値。
"""
function divergence_vector(point)
    validate_point(point)
    sum(exp, point)
end

"""
    gradient_divergence_vector(point)

ベクトル場の発散の解析的な勾配を評価する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

評価前に `validate_point` で有限な実数座標か確認する。不正な入力は `ArgumentError`。

# 返り値

x・y・z方向の成分を持つ3要素のタプル。
"""
function gradient_divergence_vector(point)
    validate_point(point)
    exp.(point)
end

"""
    laplacian_vector(point)

ベクトル場の各成分の解析的なラプラシアンを評価する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

評価前に `validate_point` で有限な実数座標か確認する。不正な入力は `ArgumentError`。

# 返り値

各成分のラプラシアンを持つ3要素のタプル。
"""
function laplacian_vector(point)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_point(point)
    x, y, z = point
    (3sin(y) * exp(2z) + exp(x), 3sin(z) * exp(2x) + exp(y), 3sin(x) * exp(2y) + exp(z))
end

end
