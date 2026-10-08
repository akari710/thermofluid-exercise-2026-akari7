# F04の数値微分演習。入力検証 → 差分TODO → 収束評価 → ベクトル恒等式の残差の順に読みます。
# 3つの差分関数を実装し、F03の解析値との比較に使います。その他の関数は提供済みです。

if !isdefined(Main, :F03VectorCalculus)
    include(joinpath(@__DIR__, "F03.jl"))
end

module F04NumericalDifferentiation

using Main.F03VectorCalculus:
    scalar_field,
    vector_field,
    gradient_scalar,
    curl_vector,
    divergence_vector,
    gradient_divergence_vector,
    laplacian_vector

export forward_difference,
    backward_difference,
    centered_difference,
    convergence_study,
    centered_partial,
    curl_gradient_residual,
    divergence_curl_residual,
    product_divergence_residual,
    curl_curl_residual,
    verify_vector_identities

"""
    validate_scalar_input(x, h)

1変数の差分に使う評価点と刻み幅を検証する。

# 引数

- `x`: 有限な実数の評価点。`Bool` は不可。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。

条件を満たさない入力は `ArgumentError` で拒否する。

# 返り値

条件を満たす場合は `nothing`。
"""
function validate_scalar_input(x, h)
    x isa Real && !(x isa Bool) && isfinite(x) ||
        throw(ArgumentError("xは有限な実数にしてください"))
    h isa Real && !(h isa Bool) && isfinite(h) && h > 0 ||
        throw(ArgumentError("hは有限な正の実数にしてください"))
    nothing
end

"""
    validate_point(point)

3次元の数値微分に使う座標を検証する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。

3要素のタプルかつ有限な実数座標を求め、`Bool` は拒否する。不正な入力は `ArgumentError`。

# 返り値

条件を満たす場合は `nothing`。
"""
function validate_point(point)
    point isa Tuple && length(point) == 3 ||
        throw(ArgumentError("点は3要素のタプルで指定してください"))
    all(value -> value isa Real && !(value isa Bool) && isfinite(value), point) ||
        throw(ArgumentError("点の座標は有限な実数にしてください"))
    nothing
end

"""
    forward_difference(f, x, h)

前進差分で1変数関数の微分を近似する。

# 引数

- `f`: 実数の評価点で値を返す関数。
- `x`: 有限な実数の評価点。`Bool` は不可。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。

不正な評価点や刻み幅は `ArgumentError` で拒否する。

# 返り値

実装後は指定点の差分商。

# 受講生のToDo

前進差分商を返す処理を実装する。配布時は `未実装 F04: forward_difference` で停止する。
"""
function forward_difference(f, x, h)
    validate_scalar_input(x, h)
    # TODO(F04): 前進差分商を実装する。
    error("未実装 F04: forward_difference")
end

"""
    backward_difference(f, x, h)

後退差分で1変数関数の微分を近似する。

# 引数

- `f`: 実数の評価点で値を返す関数。
- `x`: 有限な実数の評価点。`Bool` は不可。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。

不正な評価点や刻み幅は `ArgumentError` で拒否する。

# 返り値

実装後は指定点の差分商。

# 受講生のToDo

後退差分商を返す処理を実装する。配布時は `未実装 F04: backward_difference` で停止する。
"""
function backward_difference(f, x, h)
    validate_scalar_input(x, h)
    # TODO(F04): 後退差分商を実装する。
    error("未実装 F04: backward_difference")
end

"""
    centered_difference(f, x, h)

中心差分で1変数関数の微分を近似する。

# 引数

- `f`: 実数の評価点で値を返す関数。
- `x`: 有限な実数の評価点。`Bool` は不可。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。

不正な評価点や刻み幅は `ArgumentError` で拒否する。

# 返り値

実装後は指定点の差分商。

# 受講生のToDo

中心差分商を返す処理を実装する。配布時は `未実装 F04: centered_difference` で停止する。
"""
function centered_difference(f, x, h)
    validate_scalar_input(x, h)
    # TODO(F04): 中心差分商を実装する。
    error("未実装 F04: centered_difference")
end

"""
    convergence_study(f, derivative, x, spacings)

刻み幅を細かくし、3種類の差分誤差と観測次数を比較する。

# 引数

- `f`: 微分する関数。
- `derivative`: 独立した解析的な微分値を返す関数。
- `x`: 有限な実数の評価点。`Bool` は不可。
- `spacings`: 2要素以上で、有限な正の実数が狭義単調減少するベクトル。`Bool` は不可。

不正入力は `ArgumentError`。実行には3つの差分TODOの実装が必要。

# 返り値

`spacings` と、前進・後退・中心差分それぞれの `*_errors`、隣接誤差の比 `*_ratios`、比から求めた `*_orders` を持つ `NamedTuple`。
"""
function convergence_study(f, derivative, x, spacings)
    x isa Real && !(x isa Bool) && isfinite(x) ||
        throw(ArgumentError("xは有限な実数にしてください"))
    spacings isa AbstractVector && length(spacings) >= 2 ||
        throw(ArgumentError("spacingsには少なくとも2個の値が必要です"))
    all(h -> h isa Real && !(h isa Bool) && isfinite(h) && h > 0, spacings) ||
        throw(ArgumentError("spacingsは有限な正の実数にしてください"))
    all(index -> spacings[index] > spacings[index + 1], 1:(length(spacings) - 1)) ||
        throw(ArgumentError("spacingsは狭義単調減少にしてください"))

    exact = derivative(x)
    # 解析的な微分値を参照し、各刻み幅の誤差を集める。
    forward_errors = [abs(forward_difference(f, x, h) - exact) for h in spacings]
    backward_errors = [abs(backward_difference(f, x, h) - exact) for h in spacings]
    centered_errors = [abs(centered_difference(f, x, h) - exact) for h in spacings]
    ratios(errors) = errors[1:(end - 1)] ./ errors[2:end]
    forward_ratios = ratios(forward_errors)
    backward_ratios = ratios(backward_errors)
    centered_ratios = ratios(centered_errors)
    (;
        spacings,
        forward_errors,
        backward_errors,
        centered_errors,
        forward_ratios,
        backward_ratios,
        centered_ratios,
        forward_orders = log2.(forward_ratios),
        backward_orders = log2.(backward_ratios),
        centered_orders = log2.(centered_ratios),
    )
end

"""
    centered_partial(f, point, axis, h; differentiator = centered_difference)

1つの座標だけを変化させて、3次元関数の偏微分を近似する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。
- `f`: 3要素の座標タプルを受け取る関数。
- `axis`: 微分する方向の整数1・2・3。`Bool` は不可。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。
- `differentiator`: 1変数の差分商を返す関数。既定は学生が実装する `centered_difference`。

不正入力は `ArgumentError`。既定の差分関数を使う場合は中心差分TODOの実装が必要。

# 返り値

指定方向の偏微分の近似値。
"""
function centered_partial(f, point, axis, h; differentiator = centered_difference)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_point(point)
    h isa Real && !(h isa Bool) && isfinite(h) && h > 0 ||
        throw(ArgumentError("hは有限な正の実数にしてください"))
    axis isa Integer && !(axis isa Bool) && axis in 1:3 ||
        throw(ArgumentError("axisは1、2、3のいずれかにしてください"))
    slice(value) = f(ntuple(index -> index == axis ? value : point[index], 3))
    differentiator(slice, point[axis], h)
end

# The supplied analytic gradient is differentiated externally with the student's centered difference.
# The inner gradient is not student work in this residual calculation.
"""
    curl_gradient_residual(point, h)

解析的な勾配を数値微分し、勾配の回転の残差を求める。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。

座標は `validate_point` で検証する。実行には中心差分TODOの実装が必要。

# 返り値

回転の3成分の残差を持つタプル。内側の勾配にはF03の解析値を使う。
"""
function curl_gradient_residual(point, h)
    validate_point(point)
    (
        centered_partial(p -> gradient_scalar(p)[3], point, 2, h) -
        centered_partial(p -> gradient_scalar(p)[2], point, 3, h),
        centered_partial(p -> gradient_scalar(p)[1], point, 3, h) -
        centered_partial(p -> gradient_scalar(p)[3], point, 1, h),
        centered_partial(p -> gradient_scalar(p)[2], point, 1, h) -
        centered_partial(p -> gradient_scalar(p)[1], point, 2, h),
    )
end

# The supplied analytic curl is differentiated externally with the student's centered difference.
# The inner curl is not student work in this residual calculation.
"""
    divergence_curl_residual(point, h)

解析的な回転を数値微分し、回転の発散の残差を求める。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。

座標は `validate_point` で検証する。実行には中心差分TODOの実装が必要。

# 返り値

発散のスカラー残差。内側の回転にはF03の解析値を使う。
"""
function divergence_curl_residual(point, h)
    validate_point(point)
    sum(centered_partial(p -> curl_vector(p)[axis], point, axis, h) for axis in 1:3)
end

# Differentiate the product field numerically and evaluate the right-hand side analytically.
"""
    product_divergence_residual(point, h)

積の場の数値的な発散と、解析値による積の公式を比較する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。

座標は `validate_point` で検証する。実行には中心差分TODOの実装が必要。

# 返り値

数値的な左辺から解析的な右辺を引いたスカラー残差。
"""
function product_divergence_residual(point, h)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_point(point)
    lhs = sum(
        centered_partial(p -> scalar_field(p) * vector_field(p)[axis], point, axis, h)
        for axis in 1:3
    )
    rhs =
        sum(vector_field(point) .* gradient_scalar(point)) +
        scalar_field(point) * divergence_vector(point)
    lhs - rhs
end

"""
    curl_curl_residual(point, h)

解析的な回転を再び数値微分し、回転の回転の公式を比較する。

# 引数

- `point`: `(x, y, z)` に対応する3要素の座標タプル。入力は変更しない。
- `h`: 有限な正の実数の刻み幅。`Bool` は不可。

座標は `validate_point` で検証する。実行には中心差分TODOの実装が必要。

# 返り値

数値的な左辺と解析的な右辺の差を持つ3要素のタプル。
"""
function curl_curl_residual(point, h)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_point(point)
    lhs = (
        centered_partial(p -> curl_vector(p)[3], point, 2, h) -
        centered_partial(p -> curl_vector(p)[2], point, 3, h),
        centered_partial(p -> curl_vector(p)[1], point, 3, h) -
        centered_partial(p -> curl_vector(p)[3], point, 1, h),
        centered_partial(p -> curl_vector(p)[2], point, 1, h) -
        centered_partial(p -> curl_vector(p)[1], point, 2, h),
    )
    lhs .- (gradient_divergence_vector(point) .- laplacian_vector(point))
end

"""
    verify_vector_identities(n)

共通の内部領域で4つのベクトル恒等式の最大残差を求める。

# 引数

- `n`: 各方向の格子点数。5以上の奇数の整数で、`Bool` は不可。

不正な格子点数は `ArgumentError`。実行には中心差分TODOの実装が必要。

# 返り値

`curl_gradient`、`divergence_curl`、`product_divergence`、`curl_curl` を持つ `NamedTuple`。各値は同じ物理領域での最大絶対残差。
"""
function verify_vector_identities(n)
    n isa Integer && !(n isa Bool) && n >= 5 && isodd(n) ||
        throw(ArgumentError("nは5以上の奇数にしてください"))
    coordinates = range(-1.0, 1.0; length = Int(n))
    h = step(coordinates)
    curl_gradient = 0.0
    divergence_curl = 0.0
    product_divergence = 0.0
    curl_curl = 0.0
    # Compare maxima on a common physical region as the grid is refined.
    interior = filter(x -> abs(x) <= 0.75, coordinates[2:(end - 1)])
    for x in interior, y in interior, z in interior

        point = (x, y, z)
        curl_gradient = max(curl_gradient, maximum(abs, curl_gradient_residual(point, h)))
        divergence_curl = max(divergence_curl, abs(divergence_curl_residual(point, h)))
        product_divergence =
            max(product_divergence, abs(product_divergence_residual(point, h)))
        curl_curl = max(curl_curl, maximum(abs, curl_curl_residual(point, h)))
    end
    (; curl_gradient, divergence_curl, product_divergence, curl_curl)
end

end

if abspath(PROGRAM_FILE) == @__FILE__

    @doc """
        reference_function(x)

    直接実行する収束実験用の参照関数を評価する。

    # 引数

    - `x`: 実験の評価点。

    # 返り値

    指定点の関数値。
    """
    reference_function(x) = sin(x) * exp(x)

    @doc """
        reference_derivative(x)

    直接実行する収束実験用の解析的な微分値を評価する。

    # 引数

    - `x`: 実験の評価点。

    # 返り値

    独立した解析的な微分値。
    """
    reference_derivative(x) = exp(x) * (sin(x) + cos(x))
    spacings = [0.2, 0.1, 0.05, 0.025]
    println(
        F04NumericalDifferentiation.convergence_study(
            reference_function,
            reference_derivative,
            0.4,
            spacings,
        ),
    )
    println("n=9: ", F04NumericalDifferentiation.verify_vector_identities(9))
    println("n=17: ", F04NumericalDifferentiation.verify_vector_identities(17))
end
