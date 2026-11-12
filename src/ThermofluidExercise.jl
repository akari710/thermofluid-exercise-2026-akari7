# N05の共通API。周期添字 → 流束 → バッファ検証 → 刻み調整 → 既習課題の入口の順に読みます。
# TODOを実装して既存の入力・返り値の契約を保ち、N06以降でも再利用します。

module ThermofluidExercise

"""
    periodic_left_index(i, n)

1始まりの周期配列で左隣の添字を求める。

# 引数

- `i`: 1始まりの現在の添字。
- `n`: 配列の点数。

i,nはBoolを除く整数、n>=1、1<=i<=n。

# 返り値

実装後は 1からnまでの整数。

# 受講生のToDo

N02の周期左隣を共通化し、先頭から左へ越えたときは末尾の添字を返す。1点の場合も扱う。
この関数では引数条件の検証も学生が実装する。N05の共通APIと既習N02を参照する。配布状態では未実装エラーで停止する。
"""
function periodic_left_index(i, n)
    # TODO(N05): 周期配列の左隣添字を返す。
    if !(i isa Integer) || !(n isa Integer) ||
        i isa Bool || n isa Bool || n < 1 || !(1 <= i <= n)
        throw(ArgumentError("添字は1 <= i <= n、点数はn >= 1の整数です（Bool不可）"))
    end

    return i == 1 ? n : i - 1
    # 条件 ? A : B は，条件が真ならA、偽ならBを返す．
end

"""
    periodic_right_index(i, n)

1始まりの周期配列で右隣の添字を求める。

# 引数

- `i`: 1始まりの現在の添字。
- `n`: 配列の点数。

i,nはBoolを除く整数、n>=1、1<=i<=n。

# 返り値

実装後は 1からnまでの整数。

# 受講生のToDo

周期右隣の添字を求め、末尾から右へ越えたときは先頭へ戻す。1点の場合も扱う。
左隣APIと同じ引数条件の検証も学生が実装する。N04の周期左右隣接とN05の共通APIを参照する。配布状態では未実装エラーで停止する。
"""
function periodic_right_index(i, n)
    # TODO(N05): 周期配列の右隣添字を返す。
    if !(i isa Integer) || !(n isa Integer) ||
        i isa Bool || n isa Bool || n < 1 || !(1 <= i <= n)
        throw(ArgumentError("添字は1 <= i <= n、点数はn >= 1の整数です（Bool不可）"))
    end

    return i == n ? 1 : i + 1
end

"""
    linear_flux(u, speed)

有限な実数の線形流束を評価する。

# 引数

- `u`: 有限な実数の場の値。
- `speed`: 有限な実数の移流速度。

u,speedの有限実数を検証する。

# 返り値

実装後は スカラーの線形流束。

# 受講生のToDo

1点の場と速度を有限実数として検証し、N01で使った線形流束をスカラーで返す。配列の更新は呼出し側の責務である。
入力検証も学生の実装範囲。N05の共通APIを参照する。配布状態では未実装エラーで停止する。
"""
function linear_flux(u, speed)
    # TODO(N05): 有限実数の入力を検証し、線形流束を返す。
    if !(u isa Real && isfinite(u))
        throw(ArgumentError("uは有限な実数にしてください"))
    end
    if !(speed isa Real && isfinite(speed))
        throw(ArgumentError("speedは有限な実数にしてください"))
    end

    return speed * u
end

"""
    burgers_flux(u)

有限な実数のBurgers流束を評価する。

# 引数

- `u`: 有限な実数の場の値。

負値も受け付け、有限実数を検証する。

# 返り値

実装後は スカラーのBurgers流束。

# 受講生のToDo

1点の場を有限実数として検証し、既習N02のBurgers流束をスカラーで返す。負値も流束の入力として扱う。
入力検証も学生の実装範囲。配布状態では未実装エラーで停止する。
"""
function burgers_flux(u)
    # TODO(N05): 有限実数の入力を検証し、Burgers流束を返す。
    if !(u isa Real && isfinite(u))
        throw(ArgumentError("流束の入力を有限な実数にしてください"))
    end

    return u^2 / 2
end

"""
    validate_buffers(u_new, u_old)

新旧バッファが安全に更新できる入力か検証する。

# 引数

- `u_new`: 更新結果を書き込む配列。旧配列と重ならない独立したバッファ。
- `u_old`: 更新前の値を読む配列。更新中は変更しない。

各軸3点以上・1始まりの同形状の浮動小数配列。非aliasとは同じ記憶領域を共有しないこと。
この条件と有限な旧値を検証する。不正入力は実装後の検証でArgumentError。

# 返り値

実装後は 成功時はnothing。

# 受講生のToDo

新旧配列の型・形状・軸・点数・有限旧値と、同じ記憶領域を共有しないことを検査する。出力の初期値の有限性は要求しない。
全点で同じ旧場を読むため、入力を書き換えず、不正なら更新前にArgumentErrorを投げる検証APIを実装する。
N01〜N04のバッファ検証とN05の共通APIを参照する。配布状態では未実装エラーで停止する。
"""
function validate_buffers(u_new, u_old)
    # TODO(N05): 新旧配列の形状・型・非alias・旧値の有限性を検証する。
    if !(u_new isa AbstractArray{<:AbstractFloat} &&
        u_old isa AbstractArray{<:AbstractFloat})
        throw(ArgumentError("新旧バッファは浮動小数の配列にしてください"))
    end

    if size(u_new) != size(u_old)
        throw(ArgumentError("新旧バッファは同じ形状にしてください"))
    end

    if ndims(u_old) == 0 || any(n -> n < 3, size(u_old))
        throw(ArgumentError("新旧バッファは各軸に3点以上必要です"))
    end

    if Base.has_offset_axes(u_new, u_old)
        throw(ArgumentError("新旧バッファの添字は各軸で1始まりにしてください"))
    end

    if Base.mightalias(u_new, u_old)
        throw(ArgumentError("新旧バッファは独立した記憶領域を使ってください"))
    end

    if !all(isfinite, u_old)
        throw(ArgumentError("旧バッファの値はすべて有限にしてください"))
    end

    return nothing
end

"""
    fit_timestep(dt_max, t_final)

刻み上限を守り、最終時刻へ到達する刻みとステップ数を決める。

# 引数

- `dt_max`: 有限な正の刻み上限。
- `t_final`: 計算の最終時刻。

dt_max,t_finalは有限な正値。表現不能な刻みやステップ数も拒否する。

# 返り値

実装後は `(;steps, dt)` のNamedTuple。stepsは正整数、dtは有限正値でdt_max以下、steps回でt_finalへ到達する。

# 受講生のToDo

最終時刻へ等間隔で到達し、刻みがdt_max以下になる正整数stepsと正のdtを返す。
有限正の入力と、ステップ数・刻みの表現可能性の検証も実装する。N01〜N04の時間刻み調整とN05の共通APIを参照する。配布状態では未実装エラーで停止する。
"""
function fit_timestep(dt_max, t_final)
    # TODO(N05): 安定上限を超えない刻みとステップ数を最終時刻から返す。
    if !(dt_max isa Real && isfinite(dt_max) && dt_max > 0)
        throw(ArgumentError("dt_maxは有限な正の実数にしてください"))
    end
    if !(t_final isa Real && isfinite(t_final) && t_final > 0)
        throw(ArgumentError("t_finalは有限な正の実数にしてください"))
    end

    count = t_final / dt_max
    if !isfinite(count) || count >= typemax(Int)
        throw(ArgumentError("ステップ数が表現範囲を超えています"))
    end

    steps = max(1, ceil(Int, count))
    dt = Float64(t_final / steps)

    # 丸めによって刻みが上限を超えたら、分割数を1つ増やす。
    if dt > dt_max && steps < typemax(Int)
        steps += 1
        dt = Float64(t_final / steps)
    end

    if !(isfinite(dt) && 0 < dt <= dt_max)
        throw(ArgumentError("上限を守る正の刻みをFloat64で表現できません"))
    end

    return (; steps, dt)
end

module N01

"""
    upwind_step!(args...; kwargs...)

N01の風上更新を共通APIへ移す。

# 引数

- `args`: 新旧ベクトル、有限正の移流速度、時間刻み、格子幅の順。新配列に書込み、旧配列を保持する。
- `kwargs`: キーワード引数はない。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は更新した新配列。内部点を更新し、境界値は呼出し側で適用する。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function upwind_step!(args...; kwargs...)
    # TODO(N05): 既習のN01.upwind_step!を共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 5 && isempty(kwargs) ||
        throw(ArgumentError("新旧配列、c、dt、dxの5引数を指定してください"))
    u_new, u_old, c, dt, dx = args

    all(u -> u isa AbstractVector{<:AbstractFloat}, (u_new, u_old)) ||
        throw(ArgumentError("N01の新旧バッファは浮動小数ベクトルにしてください"))
    parentmodule(@__MODULE__).validate_buffers(u_new, u_old)

    # N01固有の入力チェック
    all(value -> value isa Real, (c, dt, dx)) ||
        throw(ArgumentError("c、dt、dxは実数にしてください"))
    all(isfinite, (c, dt, dx)) || throw(ArgumentError("c、dt、dxを有限値にしてください"))
    c > 0 || throw(ArgumentError("N01では正の移流速度だけを扱います"))
    dt > 0 || throw(ArgumentError("dtは正にしてください"))
    dx > 0 || throw(ArgumentError("dxは正にしてください"))

    courant = c * dt / dx
    copyto!(u_new, u_old)
    for i in 2:(length(u_old)-1)
        u_new[i] = u_old[i] - courant * (u_old[i] - u_old[i-1])
    end
    return u_new
end

"""
    centered_step!(args...; kwargs...)

N01の中心差分更新を共通APIへ移す。

# 引数

- `args`: 新旧ベクトル、有限正の移流速度、時間刻み、格子幅の順。新配列に書込み、旧配列を保持する。
- `kwargs`: キーワード引数はない。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は更新した新配列。内部点を更新し、境界値は呼出し側で適用する。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function centered_step!(args...; kwargs...)
    # TODO(N05): 既習のN01.centered_step!を共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 5 && isempty(kwargs) ||
        throw(ArgumentError("新旧配列、c、dt、dxの5引数を指定してください"))
    u_new, u_old, c, dt, dx = args

    all(u -> u isa AbstractVector{<:AbstractFloat}, (u_new, u_old)) ||
        throw(ArgumentError("N01の新旧バッファは浮動小数ベクトルにしてください"))
    parentmodule(@__MODULE__).validate_buffers(u_new, u_old)

    all(value -> value isa Real, (c, dt, dx)) ||
        throw(ArgumentError("c、dt、dxは実数にしてください"))
    all(isfinite, (c, dt, dx)) || throw(ArgumentError("c、dt、dxを有限値にしてください"))
    c > 0 || throw(ArgumentError("N01では正の移流速度だけを扱います"))
    dt > 0 || throw(ArgumentError("dtは正にしてください"))
    dx > 0 || throw(ArgumentError("dxは正にしてください"))

    courant = c * dt / dx
    copyto!(u_new, u_old)

    for i in 2:(length(u_old)-1)
        u_new[i] = u_old[i] - courant * (u_old[i+1] - u_old[i-1]) / 2
    end

    return u_new
end

"""
    apply_boundary!(args...; kwargs...)

N01の左流入・右流出境界を共通APIへ移す。

# 引数

- `args`: 端点を書き換えるベクトル。
- `kwargs`: `left_value` は有限な左境界値（既定1.0）。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は境界を適用した入力ベクトル。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function apply_boundary!(args...; kwargs...)
    # TODO(N05): 既習のN01.apply_boundary!を共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 1 ||
        throw(ArgumentError("境界を適用するベクトルを1つ指定してください"))
    all(key -> key === :left_value, keys(kwargs)) ||
        throw(ArgumentError("指定できるキーワード引数はleft_valueだけです"))

    u = args[1]
    left_value = get(kwargs, :left_value, 1.0)

    u isa AbstractVector{<:Real} && length(u) >= 2 ||
        throw(ArgumentError("uは2点以上の実数ベクトルにしてください"))
    left_value isa Real && isfinite(left_value) ||
        throw(ArgumentError("left_valueは有限な実数にしてください"))

    u[1] = left_value
    u[end] = u[end-1]
    return u
end
end
module N02

"""
    nonlinear_upwind_step!(args...; kwargs...)

N02の非線形風上更新を共通APIへ移す。

# 引数

- `args`: 新旧ベクトル、時間刻み、格子幅の順。新旧は独立し、旧場を保持する。
- `kwargs`: `boundary` は `:fixed` または `:periodic`（既定 `:fixed`）。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は境界も適用した新配列。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function nonlinear_upwind_step!(args...; kwargs...)
    # TODO(N05): 既習のN02.nonlinear_upwind_step!を共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 4 ||
        throw(ArgumentError("新旧配列、dt、dxの4引数を指定してください"))
    all(key -> key === :boundary, keys(kwargs)) ||
        throw(ArgumentError("指定できるキーワード引数はboundaryだけです"))

    u_new, u_old, dt, dx = args
    boundary = get(kwargs, :boundary, :fixed)
    common = parentmodule(@__MODULE__)

    # 書き込む前に、N02固有の条件と共通の配列条件を検証する。
    boundary in (:fixed, :periodic) ||
        throw(ArgumentError("boundaryは:fixedまたは:periodicです"))
    all(u -> u isa AbstractVector{<:AbstractFloat}, (u_new, u_old)) ||
        throw(ArgumentError("新旧バッファは浮動小数の一次元配列です"))
    common.validate_buffers(u_new, u_old)
    all(v -> v isa Real && isfinite(v) && v > 0, (dt, dx)) ||
        throw(ArgumentError("dtとdxは有限な正の実数です"))

    copyto!(u_new, u_old)
    n = length(u_old)
    periodic = boundary === :periodic
    indices = periodic ? (1:n) : (2:(n-1))

    for i in indices
        left = periodic ? common.periodic_left_index(i, n) : i - 1
        u_new[i] = u_old[i] - (dt / dx) * (
            common.burgers_flux(u_old[i]) - common.burgers_flux(u_old[left])
        )
    end

    all(isfinite, u_new) ||
        error("更新後に非有限値があります。更新式を確認してください")
    return u_new
end

"""
    apply_boundary!(args...; kwargs...)

N02の固定または周期境界を共通APIへ移す。

# 引数

- `args`: 境界を適用するベクトル。周期境界では格子の全点を使う。
- `kwargs`: `boundary` は `:fixed` または `:periodic`（既定 `:fixed`）。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は境界を適用した入力ベクトル。周期境界では値を変更しない。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function apply_boundary!(args...; kwargs...)
    # TODO(N05): 既習のN02.apply_boundary!を共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 1 ||
        throw(ArgumentError("境界を適用する配列を1つ指定してください"))
    all(key -> key === :boundary, keys(kwargs)) ||
        throw(ArgumentError("指定できるキーワード引数はboundaryだけです"))

    u = args[1]
    boundary = get(kwargs, :boundary, :fixed)

    boundary in (:fixed, :periodic) ||
        throw(ArgumentError("boundaryは:fixedまたは:periodicです"))
    length(u) >= 3 ||
        throw(ArgumentError("境界配列には3点以上が必要です"))

    if boundary === :fixed
        u[1] = 1.0
        u[end] = u[end-1]
    end
    return u
end
end
module N03

"""
    diffusion_step!(args...; kwargs...)

N03の熱拡散更新を共通APIへ移す。

# 引数

- `args`: 新旧ベクトル、時間刻み、格子幅、有限正の拡散係数の順。新旧は独立し、旧場を保持する。
- `kwargs`: `boundary` は固定温度の `:fixed`（既定）または断熱の `:insulated`。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は選択した境界も適用した新配列。安定上限を超える刻みも実験用に受け付ける。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function diffusion_step!(args...; kwargs...)
    # TODO(N05): 既習のN03.diffusion_step!を共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 5 ||
        throw(ArgumentError("新旧配列、dt、dx、diffusivityの5引数を指定してください"))
    all(key -> key === :boundary, keys(kwargs)) ||
        throw(ArgumentError("指定できるキーワード引数はboundaryだけです"))

    u_new, u_old, dt, dx, diffusivity = args
    boundary = get(kwargs, :boundary, :fixed)

    # 書き込む前に境界条件と新旧配列を検証する。
    boundary === :fixed || boundary === :insulated ||
        throw(ArgumentError("boundaryは:fixedまたは:insulatedです"))
    all(u -> u isa AbstractVector{<:AbstractFloat}, (u_new, u_old)) ||
        throw(ArgumentError("新旧バッファは浮動小数の一次元配列です"))
    parentmodule(@__MODULE__).validate_buffers(u_new, u_old)

    for (value, name) in ((dt, "dt"), (dx, "dx"), (diffusivity, "diffusivity"))
        value isa Real && isfinite(value) && value > 0 ||
            throw(ArgumentError("$name は有限な正の実数です"))
    end

    r = diffusivity * dt / dx^2
    isfinite(r) && r > 0 ||
        throw(ArgumentError("実効Fourier数は有限な正の値です"))

    for i in 2:(length(u_old)-1)
        u_new[i] = u_old[i] + r * (u_old[i+1] - 2 * u_old[i] + u_old[i-1])
    end

    apply_boundary!(u_new, u_old, r; boundary)
    return u_new
end

"""
    apply_boundary!(args...; kwargs...)

N03の固定温度または断熱端点の更新を共通APIへ移す。

# 引数

- `args`: 新旧ベクトル、有限正の実効Fourier数の順。新旧は独立し、旧場を保持する。
- `kwargs`: `boundary` は固定温度の `:fixed`（既定）または断熱の `:insulated`。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は両端を更新した新ベクトル。内部点は変更しない。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function apply_boundary!(args...; kwargs...)
    # TODO(N05): 既習のN03.apply_boundary!を共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 3 ||
        throw(ArgumentError("新旧配列とrの3引数を指定してください"))
    all(key -> key === :boundary, keys(kwargs)) ||
        throw(ArgumentError("指定できるキーワード引数はboundaryだけです"))

    u_new, u_old, r = args
    boundary = get(kwargs, :boundary, :fixed)

    # 両端を書き換える前に入力を検証する。
    boundary === :fixed || boundary === :insulated ||
        throw(ArgumentError("boundaryは:fixedまたは:insulatedです"))
    all(u -> u isa AbstractVector{<:AbstractFloat}, (u_new, u_old)) ||
        throw(ArgumentError("新旧バッファは浮動小数の一次元配列です"))
    parentmodule(@__MODULE__).validate_buffers(u_new, u_old)
    r isa Real && isfinite(r) && r > 0 ||
        throw(ArgumentError("rは有限な正の実数です"))

    if boundary === :fixed
        u_new[1] = 0.0
        u_new[end] = 0.0
    else
        u_new[1] = u_old[1] + 2 * r * (u_old[2] - u_old[1])
        u_new[end] = u_old[end] + 2 * r * (u_old[end-1] - u_old[end])
    end

    return u_new
end
end
module N04

"""
    stable_timestep(args...; kwargs...)

N04の合成安定条件による刻み選択を共通APIへ移す。

# 引数

- `args`: 有限非負の最大速度、有限正の格子幅、有限非負の拡散係数の順。速度と拡散係数を同時に0にしない。
- `kwargs`: `safety` は0より大きく1以下（既定0.8）。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は有限正の時間刻み。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function stable_timestep(args...; kwargs...)
    # TODO(N05): 既習のN04.stable_timestepを共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 3 ||
        throw(ArgumentError("max_speed、dx、diffusivityの3引数を指定してください"))
    all(key -> key === :safety, keys(kwargs)) ||
        throw(ArgumentError("指定できるキーワード引数はsafetyだけです"))

    max_speed, dx, diffusivity = args
    safety = get(kwargs, :safety, 0.8)

    max_speed isa Real && isfinite(max_speed) && max_speed >= 0 ||
        throw(ArgumentError("max_speedは有限な非負の実数にしてください"))
    dx isa Real && isfinite(dx) && dx > 0 ||
        throw(ArgumentError("dxは有限な正の実数にしてください"))
    diffusivity isa Real && isfinite(diffusivity) && diffusivity >= 0 ||
        throw(ArgumentError("diffusivityは有限な非負の実数にしてください"))
    safety isa Real && isfinite(safety) && 0 < safety <= 1 ||
        throw(ArgumentError("safetyは0より大きく1以下の有限な実数にしてください"))
    max_speed > 0 || diffusivity > 0 ||
        throw(ArgumentError("max_speedとdiffusivityを同時に0にはできません"))

    denominator = max_speed / dx + 2 * diffusivity / dx^2
    isfinite(denominator) && denominator > 0 ||
        throw(ArgumentError("有限な正の時間刻みを計算できません"))

    dt = safety / denominator
    isfinite(dt) && dt > 0 ||
        throw(ArgumentError("有限な正の時間刻みを計算できません"))

    return dt
end

"""
    advection_diffusion_step!(args...; kwargs...)

N04の移流拡散更新を共通APIへ移す。

# 引数

- `args`: 新旧ベクトル、時間刻み、格子幅、有限非負の拡散係数の順。新旧は独立し、旧場を保持する。
- `kwargs`: `model` を明示する。`advection=true`、`speed=1.0` は既習入口と同じ。

対応する既習関数と同じ入力検証を保つ。未実装時は未実装エラーで停止する。

# 返り値

実装後は全周期格子点を更新した新配列。

# 受講生のToDo

既習の同名関数を共通APIを使う処理へ移す。入力検証・更新範囲・入力保持を保つ。
"""
function advection_diffusion_step!(args...; kwargs...)
    # TODO(N05): 既習のN04.advection_diffusion_step!を共通APIを使う処理へ移し、既存の契約を保つ。
    length(args) == 5 ||
        throw(ArgumentError("新旧配列、dt、dx、diffusivityの5引数を指定してください"))
    all(key -> key in (:model, :advection, :speed), keys(kwargs)) ||
        throw(ArgumentError("指定できるキーワード引数はmodel、advection、speedです"))
    haskey(kwargs, :model) ||
        throw(ArgumentError("modelを明示してください"))

    u_new, u_old, dt, dx, diffusivity = args
    model = kwargs[:model]
    advection = get(kwargs, :advection, true)
    speed = get(kwargs, :speed, 1.0)
    common = parentmodule(@__MODULE__)

    model in (:linear, :nonlinear) ||
        throw(ArgumentError("modelは:linearまたは:nonlinearです"))
    advection isa Bool ||
        throw(ArgumentError("advectionはBoolにしてください"))
    all(u -> u isa AbstractVector{<:AbstractFloat}, (u_new, u_old)) ||
        throw(ArgumentError("新旧配列は浮動小数の一次元配列にしてください"))
    common.validate_buffers(u_new, u_old)

    all(value -> value isa Real && isfinite(value), (dt, dx, diffusivity, speed)) ||
        throw(ArgumentError("dt、dx、diffusivity、speedは有限な実数にしてください"))
    dt > 0 || throw(ArgumentError("dtは正にしてください"))
    dx > 0 || throw(ArgumentError("dxは正にしてください"))
    diffusivity >= 0 ||
        throw(ArgumentError("diffusivityは非負にしてください"))
    model === :linear && speed <= 0 &&
        throw(ArgumentError("線形モデルのspeedは正にしてください"))
    model === :nonlinear && advection && any(value -> value < 0, u_old) &&
        throw(ArgumentError("非線形移流では旧配列の値を非負にしてください"))

    max_speed = if !advection
        0.0
    elseif model === :linear
        speed
    else
        maximum(u_old)
    end
    cfl = max_speed * (dt / dx)
    fo = diffusivity * (dt / dx) / dx
    all(isfinite, (cfl, fo)) ||
        throw(ArgumentError("CFL・Foが有限ではありません"))
    cfl + 2 * fo <= 1 + 32 * eps(Float64) ||
        throw(ArgumentError("合成安定条件CFL+2Fo<=1を超えています"))

    n = length(u_old)
    for i in 1:n
        left = common.periodic_left_index(i, n)
        right = common.periodic_right_index(i, n)

        advection_term = if advection
            flux_i = model === :linear ?
                     common.linear_flux(u_old[i], speed) :
                     common.burgers_flux(u_old[i])
            flux_left = model === :linear ?
                        common.linear_flux(u_old[left], speed) :
                        common.burgers_flux(u_old[left])
            -(dt / dx) * (flux_i - flux_left)
        else
            0.0
        end

        diffusion_term = fo * (u_old[right] - 2 * u_old[i] + u_old[left])

        value = u_old[i] + advection_term + diffusion_term
        isfinite(value) ||
            error("更新後に非有限値があります。入力値と更新式を確認してください")
        u_new[i] = value
    end

    return u_new
end
end
include("N06Advection.jl")
include("N07Transport.jl")
include("N08N09Elliptic.jl")
end
