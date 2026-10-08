# N08/N09: 発展のSOR TODO → Gauss–Seidel入口 → 提供の反復処理。

module RelaxationExtension
using ThermofluidExercise
const E=ThermofluidExercise.Elliptic
export gauss_seidel_step!, sor_step!, solve_relaxation

"""
    sor_step!(u, f, dx, dy; omega = 1.5)

逐次更新と緩和係数でSORの1ステップを計算する。

# 引数

- `u`: 有限な浮動小数の節点場。x,y順。内部点をその場で更新する。
- `f`: Poisson方程式の右辺行列。入力は変更しない。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `omega`: 有限な実数で0より大きく2より小さい緩和係数。

u,fは同形状・1始まり・各軸3点以上で、記憶領域を共有しない。fは有限。不正入力は書込み前にArgumentError。

# 返り値

実装後は書き換えたu。fとDirichlet境界は保持する。

# 受講生のToDo

uの内部点を教材の走査順で逐次更新し、同じsweepで既に更新した隣接値も使う。四辺は書き換えず、fを保持する。
omega=1は提供のGauss–Seidel入口でも使う。入力・係数・omegaの検証は提供済み。
[N09課題「Gauss–Seidel・SORの更新条件」](https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/N09.html#relaxation-formulation)を参照する。配布状態では検証後に未実装エラーで停止する。
"""
function sor_step!(u, f, dx, dy; omega = 1.5)
    E.buffers(u, f)
    E.require(all(isfinite, u), "更新前の場は有限値が必要です")
    ax, ay, d = E.coefficients(dx, dy)
    E.require(omega isa Real && isfinite(omega) && 0 < omega < 2, "omegaは0<omega<2です")
    # TODO_BEGIN sor
    # TODO(N08-N09): 逐次更新と緩和係数を用いてSORの1ステップを実装する。
    error("未実装 N08-N09: RelaxationExtension.sor_step! (発展)")
    # TODO_END sor
    u
end

"""
    gauss_seidel_step!(u, f, dx, dy)

緩和係数1のSOR入口としてGauss–Seidel更新を行う。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `f`: Poisson方程式の右辺行列。入力は変更しない。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。

# 返り値

更新したu。SOR TODOの実装が必要。
"""
gauss_seidel_step!(u, f, dx, dy) = sor_step!(u, f, dx, dy; omega = 1.)

"""
    solve_relaxation(
        u0,
        g,
        f,
        dx,
        dy;
        method = :sor,
        omega = 1.5,
        atol = 1e-10,
        rtol = 1e-12,
        maxiter = 200000,
    )

提供の反復処理から作業用配列を逐次更新する。

# 引数

- `u0`: 初期場の行列。作業用コピーを使い、入力は変更しない。
- `g`: Dirichlet境界値の行列。入力は変更しない。
- `f`: Poisson方程式の右辺行列。入力は変更しない。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `method`: 反復法の識別子。
- `omega`: 更新の緩和係数。許容範囲は下記の検証に従う。
- `atol`: 有限で非負の絶対許容値。
- `rtol`: 有限で非負の相対許容値。
- `maxiter`: Boolを除く正整数の反復上限。

# 返り値

Elliptic.solve_driverと同じ収束情報のNamedTuple。u0,g,fは保持し、非収束を結果として返す。
"""
function solve_relaxation(
    u0,
    g,
    f,
    dx,
    dy;
    method = :sor,
    omega = 1.5,
    atol = 1e-10,
    rtol = 1e-12,
    maxiter = 200000,
)
    E.require(method in (:gauss_seidel, :sor), "method不正")
    E.require(omega isa Real && isfinite(omega) && 0 < omega < 2, "omega不正")
    # 提供の反復処理から、作業用配列を順次更新する。
    step! =
        (new, old, f, dx, dy)->(
            copyto!(new, old);
            sor_step!(new, f, dx, dy; omega = method == :gauss_seidel ? 1. : omega)
        )
    E.solve_driver(u0, g, f, dx, dy, step!, E.poisson_residual!; atol, rtol, maxiter)
end
end
