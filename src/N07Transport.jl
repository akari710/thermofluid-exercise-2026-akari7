# N07の更新API。境界とバッファの検証 → 面流束と温度更新 → 二成分Burgers更新の順に読みます。
# 数値TODO 5か所を実装します。提供検証は書込み前に呼び、旧場を保持します。

module N07
import ..ThermofluidExercise as Common
export thermal_stable_timestep,
    thermal_fluxes!, thermal_step!, burgers_stable_timestep, burgers_step!
# 提供: 形状・境界・書込み前検証。数値処理はTODOで実装する。
"""
    require(ok, msg)

条件が成立するか検証する。

# 引数

- `ok`: 確認する真偽値。
- `msg`: 条件が偽の場合のエラー文。

# 返り値

真ならnothing、偽ならArgumentError。
"""
require(ok, msg) = ok ? nothing : throw(ArgumentError(msg))

"""
    finite(v)

有限な実数か判定する。

# 引数

- `v`: 検証する実数。

# 返り値

Bool。
"""
finite(v) = v isa Real && isfinite(v)

"""
    positive(v, name)

有限な正値か検証する。

# 引数

- `v`: 検証する実数。
- `name`: 検証エラーに表示する項目名。

# 返り値

成功時はnothing。
"""
positive(v, name) = require(finite(v) && v > 0, "$name は有限正値です")
const SIDES = (:west, :east, :south, :north)

"""
    boundaries(kind::Symbol; value = 0.0)

4辺に同じ種類の境界を設定する。

# 引数

- `kind`: 境界の種類。
- `value`: 境界に設定する値。

# 返り値

検証済みの境界NamedTuple。
"""
function boundaries(kind::Symbol; value = 0.0)
    b = (; kind, value = Float64(value))
    bc = (; west = b, east = b, south = b, north = b)

    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_boundary(bc, 0., 0.)
    bc
end

"""
    channel_boundaries(wall::Symbol)

左流入・右流出と指定した上下壁の境界を作る。

# 引数

- `wall`: 上下壁の境界（:dirichletまたは:insulated）。

# 返り値

west,east,south,northの境界NamedTuple。
"""
function channel_boundaries(wall::Symbol)
    require(wall in (:dirichlet, :insulated), "上下壁は固定温度または断熱です")
    (;
        west = (; kind = :inflow, value = 1.),
        east = (; kind = :outflow, value = 0.),
        south = (; kind = wall, value = 0.),
        north = (; kind = wall, value = 0.),
    )
end

"""
    validate_boundary(bc, cx, cy)

境界の型・値・対向周期辺・速度との組合せを検証する。

# 引数

- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとFloat64のvalueを持つ。
- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。

開境界は左流入・右流出で `cx > 0, cy = 0`、閉じた非周期系は速度0。
周期辺は対向辺で指定する。不正入力は `ArgumentError`。

# 返り値

nothing。
"""
function validate_boundary(bc, cx, cy)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(bc isa NamedTuple && keys(bc) == SIDES, "bcはwest,east,south,northです")
    for b in values(bc)
        require(
            b isa NamedTuple &&
                keys(b) == (:kind, :value) &&
                b.kind in (:periodic, :dirichlet, :insulated, :inflow, :outflow) &&
                b.value isa Float64 &&
                isfinite(b.value),
            "境界の種類または値が不正です",
        )
        require(
            b.kind in (:dirichlet, :inflow) || b.value == 0.,
            "値が不要な境界はvalue=0.0です",
        )
    end
    for (a, b) in ((bc.west, bc.east), (bc.south, bc.north))
        require(
            (a.kind == :periodic) == (b.kind == :periodic),
            "周期辺は対向辺で指定します",
        )
    end
    kinds = map(b->b.kind, values(bc))
    if :inflow in kinds || :outflow in kinds
        require(
            kinds[1:2] == (:inflow, :outflow) &&
                all(k->k in (:dirichlet, :insulated), kinds[3:4]) &&
                cx > 0 &&
                cy == 0,
            "開境界は左流入・右流出、cx>0,cy=0です",
        )
    elseif !all(==(:periodic), kinds)
        require(
            all(k->k in (:periodic, :dirichlet, :insulated), kinds) && cx == cy == 0,
            "閉じた非周期系の速度は0です",
        )
    end
    nothing
end

"""
    thermal_parameters(cx, cy, kappa, dx, dy, bc; safety = 1.)

温度場の係数・格子幅・安全係数と境界を検証する。

# 引数

- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。
- `kappa`: 温度場の有限な非負拡散係数。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとFloat64のvalueを持つ。
- `safety`: 安定条件に掛ける安全係数。

# 返り値

nothing。
"""
function thermal_parameters(cx, cy, kappa, dx, dy, bc; safety = 1.)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(all(finite, (cx, cy, kappa)) && kappa >= 0, "速度は有限、kappaは有限非負です")
    positive(dx, "dx")
    positive(dy, "dy")
    positive(safety, "safety")
    require(safety <= 1, "safety<=1です")
    validate_boundary(bc, cx, cy)
end

"""
    matrix(a; old = false)

1始まりの浮動小数行列と必要な有限性を検証する。

# 引数

- `a`: 描画または検証する行列。
- `old`: 旧場として要素の有限性も検査するかを表すBool。

# 返り値

old=trueではnothing、old=falseではfalse。
"""
function matrix(a; old = false)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        a isa AbstractMatrix{<:AbstractFloat} && all(>=(3), size(a)),
        "1始まり・各軸3以上の浮動小数行列です",
    )
    require(axes(a) == map(Base.OneTo, size(a)), "1始まりの配列です")
    old && require(all(isfinite, a), "旧場は有限です")
end

"""
    independent(arrays)

複数の配列が互いに重ならないか検証する。

# 引数

- `arrays`: 互いに重ならないことを検査する配列の組。

# 返り値

nothing。
"""
function independent(arrays)
    for i in eachindex(arrays), j in (i + 1):length(arrays)
        require(!Base.mightalias(arrays[i], arrays[j]), "配列は互いに非aliasです")
    end
end

"""
    thermal_buffers(new, old)

同形状の独立した新旧温度バッファを検証する。

# 引数

- `new`: 結果を書き込む新行列。
- `old`: 読み取り専用の旧温度行列。

# 返り値

nothing。旧場は有限値を求める。
"""
function thermal_buffers(new, old)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    matrix(new)
    matrix(old; old = true)
    require(size(new) == size(old), "新旧の形状が異なります")
    independent((new, old))
end

"""
    flux_buffers(T)

x面とy面に対応する4つの流束バッファを作る。

# 引数

- `T`: 読み取り専用の旧温度行列。

# 返り値

adv_x,diff_xは(nx+1,ny)、adv_y,diff_yは(nx,ny+1)のゼロ行列を持つNamedTuple。
"""
function flux_buffers(T)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    matrix(T)
    nx, ny = size(T)
    (;
        adv_x = zeros(nx + 1, ny),
        diff_x = zeros(nx + 1, ny),
        adv_y = zeros(nx, ny + 1),
        diff_y = zeros(nx, ny + 1),
    )
end

"""
    validate_fluxes(f, T)

面流束の名前・形状・1始まりの軸と非aliasを検証する。

# 引数

- `f`: 検証する面流束のNamedTuple。
- `T`: 読み取り専用の旧温度行列。

# 返り値

nothing。温度を含む全配列が独立することを求める。
"""
function validate_fluxes(f, T)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    matrix(T; old = true)
    nx, ny = size(T)
    require(
        f isa NamedTuple && keys(f) == (:adv_x, :diff_x, :adv_y, :diff_y),
        "流束配列名が不正です",
    )
    for (a, shape) in
        zip(values(f), ((nx + 1, ny), (nx + 1, ny), (nx, ny + 1), (nx, ny + 1)))
        require(
            a isa AbstractMatrix{<:AbstractFloat} &&
                size(a) == shape &&
                axes(a) == map(Base.OneTo, shape),
            "流束配列の形状・軸が不正です",
        )
    end
    independent((T, values(f)...))
end
# 検証用の中心係数。学生はstable_timestepでこの上限の意味を説明して式を実装する。
"""
    thermal_rate(cx, cy, kappa, dx, dy, bc)

半セルの境界係数を含む合成安定条件のレートを評価する。

# 引数

- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。
- `kappa`: 温度場の有限な非負拡散係数。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとFloat64のvalueを持つ。

# 返り値

有限なスカラーレート。
"""
function thermal_rate(cx, cy, kappa, dx, dy, bc)
    ax = any(b->b.kind in (:dirichlet, :inflow), (bc.west, bc.east)) ? 3 : 2
    ay = any(b->b.kind in (:dirichlet, :inflow), (bc.south, bc.north)) ? 3 : 2
    r = abs(cx) / dx + abs(cy) / dy + kappa * (ax / dx ^ 2 + ay / dy ^ 2)

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(isfinite(r), "安定条件が表現範囲外です")
    r
end

"""
    checked_timestep(dt)

刻みをFloat64へ変換し、有限な正値か確認する。

# 引数

- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。

# 返り値

正のFloat64。
"""
function checked_timestep(dt)
    result = Float64(dt)

    # 計算や書込みの前に、入力条件をまとめて確認する。
    positive(result, "dt")
    result
end

"""
    thermal_stable_timestep(cx, cy, kappa, dx, dy, bc; safety = 0.8)

半セル境界を含む温度更新の合成上限から刻みを求める。

# 引数

- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。
- `kappa`: 温度場の有限な非負拡散係数。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとFloat64のvalueを持つ。
- `safety`: 安定条件に掛ける安全係数。

# 返り値

実装後は有限な正のFloat64。

# 受講生のToDo

TODOコメントの指示に沿って数値処理を実装する。配布状態では未実装エラーで停止する。
"""
function thermal_stable_timestep(cx, cy, kappa, dx, dy, bc; safety = 0.8)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    thermal_parameters(cx, cy, kappa, dx, dy, bc; safety)
    require(thermal_rate(cx, cy, kappa, dx, dy, bc) > 0, "全係数0では刻みを決められません")
    # TODO(N07): 半セル境界を含む合成上限から時間刻みを求める。
    dt = error("未実装 N07: thermal_stable_timestep")
    checked_timestep(dt)
end

"""
    thermal_fluxes!(fluxes, T, dx, dy; cx, cy, kappa, bc)

温度の移流・拡散流束を全ての面へ分けて書き込む。

# 引数

- `fluxes`: 書き換える面流束配列のNamedTuple。
- `T`: 読み取り専用の旧温度行列。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。
- `kappa`: 温度場の有限な非負拡散係数。
- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとFloat64のvalueを持つ。

# 返り値

実装後は書き換えたfluxes。正x・正y方向を正とし、旧温度Tは保持する。

# 受講生のToDo

TODOコメントの指示に沿って数値処理を実装する。配布状態では未実装エラーで停止する。
"""
function thermal_fluxes!(fluxes, T, dx, dy; cx, cy, kappa, bc)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    thermal_parameters(cx, cy, kappa, dx, dy, bc)
    validate_fluxes(fluxes, T)
    # TODO(N07): 正x・正y方向、移流と拡散を分けて全ての面へ書く。
    error("未実装 N07: thermal_fluxes!")
    fluxes
end

"""
    thermal_step!(Tnew, Told, dt, dx, dy; cx, cy, kappa, bc)

同じ旧温度から面流束を1回作り、温度と辺別輸送を更新する。

# 引数

- `Tnew`: 更新結果を書き込む温度行列。
- `Told`: 変更しない旧温度行列。
- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `cx`: x方向の移流速度。
- `cy`: y方向の移流速度。
- `kappa`: 温度場の有限な非負拡散係数。
- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとFloat64のvalueを持つ。

# 返り値

実装後は(;temperature=Tnew, boundary_rates=(;advective,diffusive))。辺別レートはwest,east,south,north順で流入を正とする。

# 受講生のToDo

TODOコメントの指示に沿って数値処理を実装する。配布状態では未実装エラーで停止する。
"""
function thermal_step!(Tnew, Told, dt, dx, dy; cx, cy, kappa, bc)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    thermal_parameters(cx, cy, kappa, dx, dy, bc)
    thermal_buffers(Tnew, Told)
    positive(dt, "dt")
    rate = thermal_rate(cx, cy, kappa, dx, dy, bc)
    require(dt * rate <= 1 + 32eps(Float64), "温度の安定上限を超えています")
    # TODO(N07): Common.validate_buffersも再利用し、同じ旧場の面流束を一度構築して更新と辺別レートに使う。
    error("未実装 N07: thermal_step!")
end

"""
    burgers_parameters(u, v, nu, dx, dy; safety = 1.)

二成分速度の独立性・有限性と粘性・格子幅を検証する。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `v`: 読み取り専用のy方向の旧速度行列。
- `nu`: Burgers速度の有限な非負粘性係数。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `safety`: 安定条件に掛ける安全係数。

# 返り値

成功時はnothing。
"""
function burgers_parameters(u, v, nu, dx, dy; safety = 1.)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    matrix(u; old = true)
    matrix(v; old = true)
    require(size(u) == size(v), "二成分の形状が異なります")
    independent((u, v))
    require(finite(nu) && nu >= 0, "nuは有限非負です")
    positive(dx, "dx")
    positive(dy, "dy")
    positive(safety, "safety")
    require(safety <= 1, "safety<=1です")
end

"""
    burgers_rate(u, v, nu, dx, dy)

現在の旧二成分場から局所速度の合成レートを評価する。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `v`: 読み取り専用のy方向の旧速度行列。
- `nu`: Burgers速度の有限な非負粘性係数。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。

# 返り値

有限なスカラーレート。
"""
function burgers_rate(u, v, nu, dx, dy)
    r =
        maximum(abs(u[k]) / dx + abs(v[k]) / dy for k in eachindex(u, v)) +
        2nu * (1 / dx ^ 2 + 1 / dy ^ 2)

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(isfinite(r), "安定条件が表現範囲外です")
    r
end

"""
    burgers_stable_timestep(u, v, nu, dx, dy; safety = 0.8)

旧二成分の速度と粘性から、このステップの刻み上限を求める。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `v`: 読み取り専用のy方向の旧速度行列。
- `nu`: Burgers速度の有限な非負粘性係数。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `safety`: 安定条件に掛ける安全係数。

# 返り値

実装後は正のFloat64。

# 受講生のToDo

TODOコメントの指示に沿って数値処理を実装する。配布状態では未実装エラーで停止する。
"""
function burgers_stable_timestep(u, v, nu, dx, dy; safety = 0.8)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    burgers_parameters(u, v, nu, dx, dy; safety)
    require(burgers_rate(u, v, nu, dx, dy) > 0, "全速度・粘性0では刻みを決められません")
    # TODO(N07): このステップの旧場の局所速度から合成刻みを決める。
    dt = error("未実装 N07: burgers_stable_timestep")
    checked_timestep(dt)
end

"""
    burgers_step!(unew, vnew, uold, vold, dt, dx, dy; nu)

正負の風上と拡散を使い、二成分を同じ旧場から更新する。

# 引数

- `unew`: 更新後のx速度を書き込む行列。
- `vnew`: 更新後のy速度を書き込む行列。
- `uold`: 変更しない旧x速度。
- `vold`: 変更しない旧y速度。
- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `nu`: Burgers速度の有限な非負粘性係数。

# 返り値

実装後は(;u=unew,v=vnew)。旧二成分場は変更しない。

# 受講生のToDo

TODOコメントの指示に沿って数値処理を実装する。配布状態では未実装エラーで停止する。
"""
function burgers_step!(unew, vnew, uold, vold, dt, dx, dy; nu)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    burgers_parameters(uold, vold, nu, dx, dy)
    thermal_buffers(unew, uold)
    thermal_buffers(vnew, vold)
    independent((unew, vnew, uold, vold))
    positive(dt, "dt")
    require(
        dt * burgers_rate(uold, vold, nu, dx, dy) <= 1 + 32eps(Float64),
        "Burgersの安定上限を超えています",
    )
    # TODO(N07): Common.validate_buffersと周期添字を再利用し、正負風上と拡散で二成分を同じ旧場から更新する。
    error("未実装 N07: burgers_step!")
    (; u = unew, v = vnew)
end
end
