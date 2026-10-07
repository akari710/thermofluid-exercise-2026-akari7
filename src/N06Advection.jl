# N06の2次元周期移流API。係数の検証 → 合成CFLの刻み → 旧場からの更新の順に読みます。
# stable_timestep と advection_step! の2つのTODOを実装します。

module N06
import ..ThermofluidExercise as Common
export stable_timestep, advection_step!
# 提供: 入力検証。数値式は下のTODO 2か所で実装する。
"""
    positive(v, name)

有限な正の実数を検証する。

# 引数

- `v`: 検証する値。
- `name`: 検証エラーに表示する項目名。

`v` は有限な正の実数。不正入力は `ArgumentError`。

# 返り値

成功時はtrue。
"""
function positive(v, name)
    v isa Real && isfinite(v) && v > 0 || throw(ArgumentError("$name は有限な正値です"))
end

"""
    speeds(cx, cy)

方向別の速度が有限かつ非負か検証する。

# 引数

- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。

`cx` と `cy` は有限な非負実数。不正入力は `ArgumentError`。

# 返り値

成功時はtrue。片方向の0は許可する。
"""
function speeds(cx, cy)
    all(v->v isa Real && isfinite(v) && v >= 0, (cx, cy)) ||
        throw(ArgumentError("cx,cyは有限な非負値です"))
end

"""
    stable_timestep(cx, cy, dx, dy; safety = 0.8)

2方向の合成CFLを満たす時間刻みを求める。

# 引数

- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `safety`: 安定条件に掛ける安全係数。

速度は有限非負で少なくとも一方が正、格子幅は有限正値、`0 < safety <= 1`。
条件違反や刻みの表現範囲外では `ArgumentError`。

# 返り値

実装後は有限な正のFloat64の刻み。

# 受講生のToDo

TODOコメントの指示に沿って数値処理を実装する。配布状態では未実装エラーで停止する。
"""
function stable_timestep(cx, cy, dx, dy; safety = 0.8)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    speeds(cx, cy)
    positive(dx, "dx")
    positive(dy, "dy")
    positive(safety, "safety")
    safety <= 1 && (cx > 0 || cy > 0) ||
        throw(ArgumentError("safety<=1、少なくとも一方の速度は正です"))
    # TODO(N06): 方向別寄与を足した安定条件から刻みを求める。
    dt = error("未実装 N06: stable_timestep")

    # 計算した刻みが元の型でもFloat64でも有限な正値か確認する。
    positive(dt, "合成dt（表現範囲）")
    result = Float64(dt)
    positive(result, "Float64の合成dt")
    return result
end

"""
    advection_step!(u_new, u_old, dt, dx, dy; cx = 1.0, cy = 0.5)

2方向の周期隣接を使い、同じ旧行列から全点を更新する。

# 引数

- `u_new`: 更新結果を書き込む配列。旧配列と重ならない独立したバッファ。
- `u_old`: 更新前の値を読む配列。更新中は変更しない。
- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。

新旧は同形状・1始まり・各軸3点以上の浮動小数行列で、互いに重ならない。
旧場は有限値、速度は有限非負、刻みと格子幅は有限正値で、合成CFLが1以下であること。
条件違反は書込み前に `ArgumentError`。共通APIのバッファ検証も実装済みであること。

# 返り値

実装後は更新したu_new。u[i,j]はx[i],y[j]に対応し、旧行列を保持する。

# 受講生のToDo

TODOコメントの指示に沿って数値処理を実装する。配布状態では未実装エラーで停止する。
"""
function advection_step!(u_new, u_old, dt, dx, dy; cx = 1.0, cy = 0.5)
    u_new isa AbstractMatrix && u_old isa AbstractMatrix ||
        throw(ArgumentError("新旧は行列です"))
    Common.validate_buffers(u_new, u_old)

    # 計算や書込みの前に、入力条件をまとめて確認する。
    speeds(cx, cy)
    positive(dt, "dt")
    positive(dx, "dx")
    positive(dy, "dy")
    Cx = cx * (dt / dx)
    Cy = cy * (dt / dy)
    all(isfinite, (Cx, Cy)) && Cx + Cy <= 1 + 32eps(Float64) ||
        throw(ArgumentError("合成CFL<=1を超えています: Cx=$Cx, Cy=$Cy"))
    nx, ny = size(u_old)
    # TODO(N06): 両方向の周期隣接と流束を使い、全点を旧配列だけから更新する。
    error("未実装 N06: advection_step!")
    return u_new
end
end
