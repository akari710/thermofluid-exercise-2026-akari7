# N01の実行時入力検証と公式出力の詳細を担当する教材提供ファイルです。
# 数値計算の学習対象ではなく、受講生が読解・編集する必要はありません。

using Plots
using TOML

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "results")

"""
    validate_initial_condition_inputs(x, base, plateau, plateau_start, plateau_end)

初期分布の座標と矩形範囲を検証する。

# 引数

- `x`: x方向の座標。入力は変更しない。
- `base`: 矩形領域の外側の値。
- `plateau`: 矩形領域の内側の値。
- `plateau_start`: 矩形領域の開始座標。
- `plateau_end`: 矩形領域の終了座標。

座標は空でない有限値の狭義単調増加ベクトル、パラメータは有限実数で、矩形領域を計算領域内に置く。
不正入力は `ArgumentError`。

# 返り値

検証に成功すると `nothing`。入力は変更しない。
"""
function validate_initial_condition_inputs(x, base, plateau, plateau_start, plateau_end)
    isempty(x) && throw(ArgumentError("xを空にすることはできません"))
    all(isfinite, x) || throw(ArgumentError("xの全要素を有限値にしてください"))
    all(isfinite, (base, plateau, plateau_start, plateau_end)) ||
        throw(ArgumentError("初期条件のパラメータを有限値にしてください"))
    all(diff(x) .> 0) || throw(ArgumentError("xを狭義単調増加にしてください"))
    first(x) <= plateau_start <= plateau_end <= last(x) ||
        throw(ArgumentError("矩形領域を計算領域内に置いてください"))
    return nothing
end

"""
    validate_step_inputs(u_new, u_old, c, dt, dx)

新旧バッファの独立性と更新係数を検証する。

# 引数

- `u_new`: 更新結果を書き込む配列。旧配列と重ならない独立したバッファ。
- `u_old`: 更新前の値を読む配列。更新中は変更しない。
- `c`: 正の移流速度。
- `dt`: 時間刻み。指定可能な範囲は下記の検証に従う。
- `dx`: x方向の有限な正の格子幅。

新旧は別のオブジェクトで同長・3点以上、旧場は有限値。`c`、`dt`、`dx` は有限正値。
不正入力は書込み前に `ArgumentError`。

# 返り値

検証に成功すると `nothing`。入力は変更しない。
"""
function validate_step_inputs(u_new, u_old, c, dt, dx)
    u_new === u_old && throw(ArgumentError("新旧で別々のバッファを使ってください"))
    length(u_new) == length(u_old) >= 3 ||
        throw(ArgumentError("二つのバッファを同じ長さの3点以上にしてください"))
    all(isfinite, u_old) || throw(ArgumentError("u_oldの全要素を有限値にしてください"))
    all(isfinite, (c, dt, dx)) || throw(ArgumentError("c、dt、dxを有限値にしてください"))
    c > 0 || throw(ArgumentError("N01では正の移流速度だけを扱います"))
    dt > 0 || throw(ArgumentError("dtは正にしてください"))
    dx > 0 || throw(ArgumentError("dxは正にしてください"))
    return nothing
end

"""
    validate_boundary_inputs(u, left_value)

境界配列と左端の値を検証する。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `left_value`: 左境界へ設定する値。

配列は2点以上、左端の値は有限値。不正入力は書込み前に `ArgumentError`。

# 返り値

検証に成功すると `nothing`。入力は変更しない。
"""
function validate_boundary_inputs(u, left_value)
    length(u) >= 2 || throw(ArgumentError("uには2点以上が必要です"))
    isfinite(left_value) || throw(ArgumentError("left_valueを有限値にしてください"))
    return nothing
end

"""
    validate_simulation_inputs(scheme, nx, c, cfl, t_final)

差分法、格子点数、CFLと最終時刻を検証する。

# 引数

- `scheme`: 使用する差分法の識別子。
- `nx`: x方向の格子点数。
- `c`: 正の移流速度。
- `cfl`: 指定するCFL数。
- `t_final`: 計算の最終時刻。

`scheme` は `:upwind` または `:centered`、`nx` はBoolを除く3以上の整数。
速度・CFL・最終時刻は有限正値。不正入力は `ArgumentError`。

# 返り値

検証に成功すると `nothing`。入力は変更しない。
"""
function validate_simulation_inputs(scheme, nx, c, cfl, t_final)
    scheme in (:upwind, :centered) ||
        throw(ArgumentError("schemeには:upwindまたは:centeredを指定してください"))
    nx isa Bool && throw(ArgumentError("nxには格子点数を表す整数を指定してください"))
    nx >= 3 || throw(ArgumentError("nxは3以上にしてください"))
    all(isfinite, (c, cfl, t_final)) ||
        throw(ArgumentError("c、cfl、t_finalを有限値にしてください"))
    c > 0 || throw(ArgumentError("N01では正の移流速度だけを扱います"))
    0 < cfl <= 1 || throw(ArgumentError("cflは0 < cfl <= 1を満たす必要があります"))
    t_final > 0 || throw(ArgumentError("t_finalは正にしてください"))
    return nothing
end

"""
    write_summary(output_dir::AbstractString, upwind, centered)

2手法の診断量を機械可読なTOMLへ保存する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `upwind`: 風上差分の計算結果。
- `centered`: 中心差分の計算結果。

# 返り値

保存した `summary.toml` のパス。出力先を作成し、同名ファイルを上書きする。
"""
function write_summary(output_dir::AbstractString, upwind, centered)
    mkpath(output_dir)
    path = joinpath(output_dir, "summary.toml")

    # 計算値と条件を、保存する診断情報にまとめる。
    summary = Dict(
        "course_id" => "N01",
        "grid" => Dict("nx" => length(upwind.x), "dx" => upwind.dx),
        "upwind" => summary_section("upwind-euler", upwind),
        "centered_euler" => summary_section("centered-euler", centered),
    )
    open(path, "w") do io
        TOML.print(io, summary; sorted = true)
    end
    return path
end

"""
    make_plots(output_dir::AbstractString, upwind, centered)

初期値と最終値を比較する2つの公式図を作る。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `upwind`: 風上差分の計算結果。
- `centered`: 中心差分の計算結果。

# 返り値

`upwind` と `centered` のPNGパスを持つ `NamedTuple`。出力先を作成し図を保存する。
"""
function make_plots(output_dir::AbstractString, upwind, centered)
    mkpath(output_dir)
    upwind_path = joinpath(output_dir, "upwind.png")
    centered_path = joinpath(output_dir, "centered-euler.png")

    for (result, title, path) in (
        (upwind, "Upwind + Euler (stable)", upwind_path),
        (centered, "Centered + Euler (intentionally unstable)", centered_path),
    )
        plot(
            result.x,
            result.u0;
            label = "Initial",
            linewidth = 2,
            xlabel = "x",
            ylabel = "u",
            title = title,
        )
        plot!(result.x, result.u; label = "Final", linewidth = 2)
        savefig(path)
    end

    # 場と診断量をまとめて返し、呼出し側で比較や保存に使えるようにする。
    return (upwind = upwind_path, centered = centered_path)
end
