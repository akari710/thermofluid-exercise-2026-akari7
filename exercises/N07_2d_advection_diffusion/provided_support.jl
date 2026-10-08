# N07: 境界条件・解析解・保存schema・来歴・容量検査と復元付き出力を提供します。

# 提供: N07セル中心格子と辺順。計算・保存・解析・作図で共有する。
using HDF5, SHA, TOML
const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "results")
const OUTPUT_NAMES = (
    "temperature.h5",
    "burgers.h5",
    "summary.toml",
    "temperature_comparison.png",
    "boundary_heat.png",
    "burgers_fields.png",
    "convergence.png",
    "plots.toml",
)
const OFFICIAL_GRIDS = ((40, 30), (80, 60), (160, 120))
const SAVE_TIMES = [0., 0.25, 0.5, 0.75, 1.]
const SIDES = ("west", "east", "south", "north")

"""
    case_id(nx, ny)

格子点数から保存ケース名を作る。

# 引数

- `nx`: x方向の格子点数。
- `ny`: y方向の格子点数。

# 返り値

n040x030形式のString。
"""
case_id(nx, ny) = "n$(lpad(nx,3,'0'))x$(lpad(ny,3,'0'))"

"""
    cell_centers(n, L)

周期終点を重複させないセル中心座標を作る。

# 引数

- `n`: 配列の点数。
- `L`: 周期領域の長さ。

# 返り値

新しい座標ベクトル。
"""
cell_centers(n, L) = collect(((1:n) .- 0.5) .* (L / n))

"""
    require(ok, msg)

保存と解析に必要な条件を検証する。

# 引数

- `ok`: 確認する真偽値。
- `msg`: 条件が偽の場合のエラー文。

# 返り値

真ならnothing、偽ならArgumentError。
"""
require(ok, msg) = ok ? nothing : throw(ArgumentError(msg))

"""
    schedule(time, t_final = 1.)

保存時刻を検証してFloat64へ変換する。

# 引数

- `time`: 保存時刻列。
- `t_final`: 計算の最終時刻。

# 返り値

新しい時刻ベクトル。
"""
function schedule(time, t_final = 1.)
    require(
        time isa AbstractVector &&
            length(time) >= 2 &&
            all(t->t isa Real && isfinite(t), time) &&
            first(time) == 0 &&
            last(time) == t_final &&
            all(diff(time) .> 0),
        "timeは0始まり・有限・狭義増加・末尾t_finalです",
    )
    Float64.(time)
end

"""
    validate_budget(time, heat, adv, diff)

熱量と辺別区間積分の形状・時刻・有限性を検証する。

# 引数

- `time`: 保存時刻列。
- `heat`: 時刻に対応する有限な熱量列。
- `adv`: west,east,south,north順の移流熱輸送積分（4,nt-1）。
- `diff`: 同じ辺順の拡散熱輸送積分（4,nt-1）。

# 返り値

nothing。
"""
function validate_budget(time, heat, adv, diff)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(time isa AbstractVector && length(time) >= 2, "timeは2点以上です")
    schedule(time, last(time))
    nt = length(time)
    require(
        heat isa AbstractVector &&
            length(heat) == nt &&
            all(v->v isa Real && isfinite(v), heat),
        "heatは時刻に対応する有限列です",
    )
    for a in (adv, diff)
        require(
            a isa AbstractMatrix &&
                size(a) == (4, nt - 1) &&
                all(v->v isa Real && isfinite(v), a),
            "熱輸送積分は有限な(4,nt-1)配列です",
        )
    end
end
const PERIODIC_CASES = ("periodic_advection", "periodic_diffusion", "periodic_combined")
const THERMAL_CASES = (
    PERIODIC_CASES...,
    "closed_fixed",
    "closed_insulated",
    "channel_fixed",
    "channel_insulated",
    ("comparison_" * split(c, "_")[2] for c in PERIODIC_CASES)...,
)

"""
    thermal_config(id)

公式温度ケースの係数・境界と初期条件を設定する。

# 引数

- `id`: 格子点数と対応するケース名。

# 返り値

cx,cy,kappa,bc,initialのNamedTuple。未知ケースは拒否する。
"""
function thermal_config(id)
    b(k, v = 0.) = (; kind = k, value = Float64(v))
    allbc(k) = (; west = b(k), east = b(k), south = b(k), north = b(k))

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(id in THERMAL_CASES, "未知の温度ケース: $id")
    if startswith(id, "periodic_") || startswith(id, "comparison_")
        mode = split(id, "_")[2]
        cx, cy = mode == "diffusion" ? (0., 0.) : (1., 0.5)
        return (;
            cx,
            cy,
            kappa = mode == "advection" ? 0. : 0.05,
            bc = allbc(:periodic),
            initial = "fourier_v1",
        )
    elseif startswith(id, "closed_")
        return (;
            cx = 0.,
            cy = 0.,
            kappa = 0.05,
            bc = allbc(id == "closed_fixed" ? :dirichlet : :insulated),
            initial = "sine_wall_v1",
        )
    end
    wall = id == "channel_fixed" ? :dirichlet : :insulated
    (;
        cx = 1.,
        cy = 0.,
        kappa = 0.05,
        bc = (;
            west = b(:inflow, 1.),
            east = b(:outflow),
            south = b(wall),
            north = b(wall),
        ),
        initial = "uniform_0.2",
    )
end

"""
    exact_temperature(x, y, t, id)

解析式を持つ温度ケースの場を評価する。

# 引数

- `x`: x方向の座標。入力は変更しない。
- `y`: y方向の座標。入力は変更しない。
- `t`: 評価する時刻。
- `id`: 格子点数と対応するケース名。

# 返り値

新しい温度行列。解析式のないケース/時刻はnothing。
"""
function exact_temperature(x, y, t, id)
    c = thermal_config(id)
    if c.initial == "fourier_v1"
        X = x .- c.cx * t
        Y = y .- c.cy * t
        k = c.kappa
        return [
            1 +
            0.2exp(-k * pi ^ 2 * t) * sin(pi * xi) +
            0.3exp(-4k * pi ^ 2 * t) * cos(2pi * yj) +
            0.1exp(-5k * pi ^ 2 * t) * sin(pi * xi + 2pi * yj) for xi in X, yj in Y
        ]
    elseif id == "closed_fixed" || (id == "closed_insulated" && t == 0)
        return [
            sin(pi * xi / 2) * sin(pi * yj) * exp(-c.kappa * ((pi / 2) ^ 2 + pi ^ 2) * t)
            for xi in x, yj in y
        ]
    elseif startswith(id, "channel_") && t == 0
        return fill(0.2, length(x), length(y))
    end
    nothing
end

"""
    exact_burgers(x, y, t; nu = 0.05)

Cole–Hopfの明示微分から二成分の解析解を求める。

# 引数

- `x`: x方向の座標。入力は変更しない。
- `y`: y方向の座標。入力は変更しない。
- `t`: 評価する時刻。
- `nu`: Burgers速度の有限な非負粘性係数。

# 返り値

(u,v)の新しい行列のタプル。解析値を数値差分で作らない。
"""
function exact_burgers(x, y, t; nu = 0.05)
    X = x .- 0.6t
    Y = y .+ 0.3t
    A = 0.2exp(-nu * pi ^ 2 * t)
    B = 0.2exp(-4nu * pi ^ 2 * t)
    phi = [1 + A * cos(pi * xi) + B * cos(2pi * yj) for xi in X, yj in Y]
    u = [
        0.6 + 2nu * A * pi * sin(pi * xi) / phi[i, j] for
        (i, xi) in enumerate(X), (j, yj) in enumerate(Y)
    ]
    v = [
        -0.3 + 4nu * B * pi * sin(2pi * yj) / phi[i, j] for
        (i, xi) in enumerate(X), (j, yj) in enumerate(Y)
    ]
    u, v
end

"""
    file_sha(path)

ファイルのSHA256を求める。

# 引数

- `path`: 対象ファイルのパス。

# 返り値

16進数String。
"""
file_sha(path) = bytes2hex(sha256(read(path)))

"""
    write_toml(path, doc)

一時出力先へ解析辞書を保存する。

# 引数

- `path`: 対象ファイルのパス。
- `doc`: 保存するTOML辞書。

# 返り値

TOML.printの結果。pathを上書きする。
"""
write_toml(path, doc) = open(io->TOML.print(io, doc; sorted = true), path, "w")

"""
    check_course_capacity(output_dir, new_total)

他課題と追加ファイルを含む全課題の容量を検査する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。
- `new_total`: 更新後に保持する成果物の容量合計。

# 返り値

成功時はnothing。標準のexercises配下で100MiBを超えると拒否する。
"""
function check_course_capacity(output_dir, new_total)
    # Find the actual student's exercises root, including custom nested result files.
    task = dirname(abspath(output_dir))
    exercises = dirname(task)
    basename(exercises) == "exercises" || return
    total = new_total
    for entry in readdir(exercises; join = true)
        entry == task && continue
        results = joinpath(entry, "results")
        isdir(results) || continue
        for (dir, _, files) in walkdir(results), name in files
            total+=filesize(joinpath(dir, name))
        end
    end
    legacy = joinpath(dirname(exercises), "results")
    if isdir(legacy)
        for (dir, _, files) in walkdir(legacy), name in files
            total+=filesize(joinpath(dir, name))
        end
    end

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(total <= 100 * 1024 ^ 2, "全課題の100 MiB上限を超えています")
end

"""
    publish(
        stage,
        output_dir,
        names;
        copy_file = (a, b)->cp(a, b; force = true),
        restore_file = (a, b)->cp(a, b; force = true),
        file_limit = 5 * 1024 ^ 2,
        task_limit = 10 * 1024 ^ 2,
    )

公式出力の容量を検査し、バックアップを作って反映する。

# 引数

- `stage`: 新しい成果物を生成する一時ディレクトリ。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `names`: 反映する公式ファイルの相対パスの列。
- `copy_file`: 新しいファイルを公式パスへ反映する関数。
- `restore_file`: バックアップから既存ファイルを復元する関数。
- `file_limit`: 1ファイルの容量上限（byte）。
- `task_limit`: 内容IDごとの容量上限（byte）。

# 返り値

nothing。反映失敗は部分書込みを含め復元する。復元失敗ではバックアップを残し場所と元エラーを示す。追加ファイルは保持する。
"""
function publish(
    stage,
    output_dir,
    names;
    copy_file = (a, b)->cp(a, b; force = true),
    restore_file = (a, b)->cp(a, b; force = true),
    file_limit = 5 * 1024 ^ 2,
    task_limit = 10 * 1024 ^ 2,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        all(n->n in OUTPUT_NAMES, names) && length(unique(names)) == length(names),
        "反映できる名前はN07公式出力だけです",
    )
    sizes = [filesize(joinpath(stage, n)) for n in names]
    require(all(sizes .<= file_limit), "出力の1ファイル上限を超えています")
    # Include retained user files in the one N07 total.
    totals = Dict{String,Int}()
    group(n) = "N07"
    for (n, size) in zip(names, sizes)
        totals[group(n)] = get(totals, group(n), 0) + size
    end
    if isdir(output_dir)
        for (dir, _, files) in walkdir(output_dir), file in files
            path = joinpath(dir, file)
            rel = relpath(path, output_dir)
            rel in names && continue
            require(
                filesize(path) <= file_limit,
                "追加ファイルの1ファイル上限を超えています",
            )
            key = group(rel)
            totals[key] = get(totals, key, 0) + filesize(path)
        end
    end
    require(all(values(totals) .<= task_limit), "出力の内容ID合計上限を超えています")
    check_course_capacity(output_dir, sum(values(totals)))

    # 新しい出力と既存ファイルのバックアップを別々に用意する。
    backup = mktempdir(; cleanup = false)
    existed = Dict{String,Bool}()
    applied = String[]
    try
        for name in names
            target = joinpath(output_dir, name)
            require(
                !ispath(target) || isfile(target),
                "出力先が通常ファイルではありません: $target",
            )
            existed[name] = isfile(target)
            if existed[name]
                saved = joinpath(backup, name)
                mkpath(dirname(saved))
                cp(target, saved)
            end
        end
        for name in names
            target = joinpath(output_dir, name)
            mkpath(dirname(target))
            push!(applied, name) # Include partial writes from a failing copy.
            copy_file(joinpath(stage, name), target)
        end
    catch original
        failures = String[]
        for name in reverse(applied)
            target = joinpath(output_dir, name)
            try
                existed[name] ? restore_file(joinpath(backup, name), target) :
                rm(target; force = true)
            catch e
                push!(failures, "$target: $e")
            end
        end
        if !isempty(failures)
            error(
                "反映と復元に失敗しました。バックアップ: $backup\n" *
                join(failures, "\n") *
                "\n元のエラー: $original",
            )
        end
        rm(backup; recursive = true)
        rethrow()
    end
    rm(backup; recursive = true)
    nothing
end

"""
    staged(action, output_dir, names; kwargs...)

一時領域で全生成・検査を成功させてから公式出力を反映する。

# 引数

- `action`: 一時ディレクトリで生成・検査を行う関数。
- `output_dir`: 公式成果物を書き出すディレクトリ。
- `names`: 反映する公式ファイルの相対パスの列。
- `kwargs`: 入口へ渡すキーワード引数。

# 返り値

nothing。生成・検査失敗では既存出力を保持する。
"""
function staged(action, output_dir, names; kwargs...)
    mktempdir() do temporary
        action(temporary)
        publish(temporary, output_dir, names; kwargs...)
    end
end

const CASE_ATTRIBUTES = (
    "nx",
    "ny",
    "dx",
    "dy",
    "grid_location",
    "cx",
    "cy",
    "kappa",
    "nu",
    "initial_condition",
    "safety",
    "timestep_policy",
    ("$(s)_$(k)" for s in SIDES for k in ("kind", "value"))...,
)
const ROOT_KEYS =
    ("schema_version", "task_id", "run_id", "equation_family", "code_sha256", "git_head")
const COMMON_DATASETS = ("x", "y", "time", "step_time", "step_dt", "save_step_index")

"""
    field_keys(family)

方程式系に応じた場のdataset名を返す。

# 引数

- `family`: temperatureまたはburgersの方程式系。

# 返り値

temperatureは1要素、burgersは(u,v)のタプル。
"""
field_keys(family) = family == "temperature" ? ("temperature",) : ("u", "v")

"""
    dataset_keys(family)

共通データ・場・熱輸送に必要なdataset名をまとめる。

# 引数

- `family`: temperatureまたはburgersの方程式系。

# 返り値

dataset名のタプル。
"""
dataset_keys(family) = (
    COMMON_DATASETS...,
    field_keys(family)...,
    (
        family == "temperature" ?
        ("advective_heat_integral", "diffusive_heat_integral") : ()
    )...,
)

"""
    source_metadata(run_id, family)

同じ実行を識別するコードhashとGit来歴をまとめる。

# 引数

- `run_id`: 同じ実行の保存場を対応付ける識別子。
- `family`: temperatureまたはburgersの方程式系。

# 返り値

HDF5ルート属性の辞書。
"""
function source_metadata(run_id, family)
    root = normpath(joinpath(@__DIR__, "..", ".."))
    paths = vcat(
        [joinpath("src", n) for n in readdir(joinpath(root, "src")) if endswith(n, ".jl")],
        [
            relpath(joinpath(@__DIR__, n), root) for
            n in ("simulate.jl", "provided_support.jl")
        ],
    )
    hashes = Dict(n => file_sha(joinpath(root, n)) for n in paths)
    code = sprint(io->TOML.print(io, hashes; sorted = true))
    head = try
        readchomp(pipeline(`git -C $root rev-parse HEAD`; stderr = devnull))
    catch
        "unversioned"
    end
    Dict{String,Any}(
        "schema_version" => 1,
        "task_id" => "N07",
        "run_id" => run_id,
        "equation_family" => family,
        "code_sha256" => code,
        "git_head" => head,
    )
end

"""
    grids_for(id)

公式ケースに必要な格子集合を選ぶ。

# 引数

- `id`: 格子点数と対応するケース名。

# 返り値

周期ケースは公式3格子、その他は80x60だけのタプル。
"""
function grids_for(id)
    id in PERIODIC_CASES || id == "periodic_burgers" ? OFFICIAL_GRIDS : ((80, 60),)
end

"""
    expected_cases(family)

方程式系ごとの公式ケースと格子の対応を作る。

# 引数

- `family`: temperatureまたはburgersの方程式系。

# 返り値

ケース名から格子名集合への辞書。
"""
function expected_cases(family)
    Dict(
        id => Set(case_id(n...) for n in grids_for(id)) for
        id in (family == "temperature" ? THERMAL_CASES : ("periodic_burgers",))
    )
end

"""
    case_metadata(id, nx, ny)

公式ケースのセル中心格子・係数・境界・刻み方針を作る。

# 引数

- `id`: 格子点数と対応するケース名。
- `nx`: x方向の格子点数。
- `ny`: y方向の格子点数。

# 返り値

ケース属性の辞書。
"""
function case_metadata(id, nx, ny)
    b = id == "periodic_burgers"
    c =
        b ?
        (;
            cx = 0.6,
            cy = -0.3,
            kappa = 0.,
            bc = thermal_config("periodic_advection").bc,
            initial = "cole_hopf_v1",
        ) : thermal_config(id)

    # 計算値と条件を、保存する診断情報にまとめる。
    doc = Dict{String,Any}(
        "nx" => nx,
        "ny" => ny,
        "dx" => 2 / nx,
        "dy" => 1 / ny,
        "grid_location" => "cell_center",
        "cx" => c.cx,
        "cy" => c.cy,
        "kappa" => c.kappa,
        "nu" => b ? 0.05 : 0.,
        "initial_condition" => c.initial,
        "safety" => 0.8,
        "timestep_policy" =>
            b ? "dynamic_old_field" :
            startswith(id, "comparison_") ? "common_combined" : "case_stable",
    )
    for side in SIDES
        v = getproperty(c.bc, Symbol(side))
        doc[side * "_kind"] = string(v.kind)
        doc[side * "_value"] = v.value
    end
    doc
end

"""
    read_attribute(o, k)

HDF5の必須属性を読む。

# 引数

- `o`: 属性を読むHDF5オブジェクト。
- `k`: 属性名。

指定した属性が必要。欠落時は `ArgumentError`。

# 返り値

属性値。欠落はArgumentError。
"""
read_attribute(o, k) =
    haskey(attributes(o), k) ? read(attributes(o)[k]) :
    throw(ArgumentError("必要属性がありません: $k"))

"""
    validate_case(c, id, grid, family)

セル中心座標・場・ステップ記録・保存時刻と熱輸送を検証する。

# 引数

- `c`: 格子・座標・条件・保存場を持つケース辞書。
- `id`: 格子点数と対応するケース名。
- `grid`: 格子点数に対応する保存グループ名。
- `family`: temperatureまたはburgersの方程式系。

# 返り値

温度ケースではnothing、Burgersではnothing。入力を変更せず、条件・形状・安定刻みの不一致を拒否する。
"""
function validate_case(c, id, grid, family)
    nx, ny = c["nx"], c["ny"]

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        nx isa Integer && ny isa Integer && nx >= 3 && ny >= 3 && grid == case_id(nx, ny),
        "$id/$grid: 格子数が不正です",
    )
    expected = case_metadata(id, nx, ny)
    require(
        all(c[k] == expected[k] for k in CASE_ATTRIBUTES),
        "$id: 境界・係数・格子配置・刻み方針が不一致です",
    )
    for (key, n, L) in (("x", nx, 2.), ("y", ny, 1.))
        require(
            c[key] isa Vector{Float64} &&
                length(c[key]) == n &&
                all(isfinite, c[key]) &&
                isapprox(c[key], cell_centers(n, L); rtol = 1e-13, atol = 1e-14),
            "$id/$key: セル中心座標が不正です",
        )
    end
    time = schedule(c["time"])
    require(time == SAVE_TIMES, "公式保存時刻が不一致です")
    nt = length(time)
    for key in field_keys(family)
        require(
            c[key] isa Array{Float64,3} &&
                size(c[key]) == (nx, ny, nt) &&
                all(isfinite, c[key]),
            "$id/$key: shape・有限値が不正です",
        )
    end
    ts, dt, idx = c["step_time"], c["step_dt"], c["save_step_index"]
    require(
        ts isa Vector{Float64} &&
            dt isa Vector{Float64} &&
            length(ts) == length(dt) > 0 &&
            all(isfinite, ts) &&
            all(v->isfinite(v) && v > 0, dt),
        "$id: step配列が不正です",
    )
    require(
        first(ts) == 0 &&
            all(diff(ts) .> 0) &&
            all(
                isapprox.(
                    ts[2:end],
                    ts[1:(end - 1)] + dt[1:(end - 1)];
                    atol = 1e-13,
                    rtol = 1e-13,
                ),
            ),
        "$id: step時刻とdtが不一致です",
    )
    require(
        idx isa AbstractVector{<:Integer} &&
            length(idx) == nt &&
            first(idx) == 0 &&
            last(idx) == length(dt) &&
            all(diff(idx) .> 0),
        "$id: 保存step対応が不正です",
    )
    require(
        all(
            isapprox(time[k], ts[idx[k]] + dt[idx[k]]; atol = 1e-13, rtol = 1e-13) for
            k in 2:nt
        ),
        "$id: 保存時刻とstepが不一致です",
    )
    if family == "temperature"
        for key in ("advective_heat_integral", "diffusive_heat_integral")
            require(
                c[key] isa Matrix{Float64} &&
                    size(c[key]) == (4, nt - 1) &&
                    all(isfinite, c[key]),
                "$id/$key: 辺・区間・有限値が不正です",
            )
        end
        # Independent metadata stability check; no numerical package is imported for reading.
        ax =
            any(c[s * "_kind"] in ("dirichlet", "inflow") for s in ("west", "east")) ? 3 : 2
        ay =
            any(c[s * "_kind"] in ("dirichlet", "inflow") for s in ("south", "north")) ? 3 :
            2
        rate =
            abs(c["cx"]) / c["dx"] +
            abs(c["cy"]) / c["dy"] +
            c["kappa"] * (ax / c["dx"] ^ 2 + ay / c["dy"] ^ 2)
        startswith(id, "comparison_") &&
            (rate = 1 / c["dx"] + 0.5 / c["dy"] + 0.1 * (1 / c["dx"] ^ 2 + 1 / c["dy"] ^ 2))
        require(all(dt .* rate .<= 0.8 + 32eps()), "$id: 安定刻み超過です")
    end
end

"""
    write_fields(path, cases, meta; official = true)

方程式系のケースをHDF5へ保存し、読取りで再検証する。

# 引数

- `path`: 対象ファイルのパス。
- `cases`: ケース名をキーにした計算結果の辞書。
- `meta`: HDF5のルート属性。
- `official`: 公式ケースと格子の集合も検査するか。

# 返り値

検証済みのmetadataとcases。Julia場は(nx,ny,nt)、属性axis_orderはtime,y,x。辺積分のJulia軸は(side,interval)、流入を正とする。
"""
function write_fields(path, cases, meta; official = true)
    family = meta["equation_family"]
    h5open(path, "w") do h
        for (k, v) in meta
            attributes(h)[k] = v
        end
        cg = create_group(h, "cases")
        for (id, grids) in sort(collect(cases); by = first)
            ig = create_group(cg, id)
            for (grid, c) in sort(collect(grids); by = first)
                g = create_group(ig, grid)
                for k in CASE_ATTRIBUTES
                    attributes(g)[k] = c[k]
                end
                for k in dataset_keys(family)
                    g[k] = c[k]
                    k != "save_step_index" && (attributes(g[k])["units"] = "1")
                    k in field_keys(family) && (attributes(g[k])["axis_order"] = "time,y,x")
                    if endswith(k, "heat_integral")
                        attributes(g[k])["axis_order"] = "interval,side"
                        attributes(g[k])["side_order"] = "west,east,south,north"
                        attributes(g[k])["positive_direction"] = "inward"
                    end
                end
            end
        end
    end
    read_fields(path; official)
end

"""
    read_fields(path; official = true)

HDF5のschema・軸・保存値・来歴と公式ケース集合を検証する。

# 引数

- `path`: 対象ファイルのパス。
- `official`: 公式ケースと格子の集合も検査するか。

# 返り値

metadataとcasesを持つNamedTuple。欠落・破損・不一致はファイル名付きエラー。
"""
function read_fields(path; official = true)
    isfile(path) || error("HDF5欠落: $(path)。simulate.jlを実行してください")
    try
        h5open(path, "r") do h
            meta = Dict(k => read_attribute(h, k) for k in ROOT_KEYS)
            family = meta["equation_family"]
            require(
                meta["schema_version"] == 1 &&
                    meta["task_id"] == "N07" &&
                    family in ("temperature", "burgers"),
                "root schema不一致",
            )
            require(
                meta["run_id"] isa String && occursin(r"^[0-9a-f]{32}$", meta["run_id"]),
                "run_id不正",
            )
            hashes = TOML.parse(meta["code_sha256"])
            require(
                !isempty(hashes) &&
                    all(
                        v->v isa String && occursin(r"^[0-9a-f]{64}$", v),
                        values(hashes),
                    ) &&
                    meta["git_head"] isa String &&
                    !isempty(meta["git_head"]),
                "コード来歴が不正です",
            )
            require(haskey(h, "cases"), "cases欠落")
            cases = Dict{String,Any}()
            for id in keys(h["cases"])
                require(id in keys(expected_cases(family)), "未知ケースです: $id")
                grids = Dict{String,Any}()
                for grid in keys(h["cases/$id"])
                    g = h["cases/$id/$grid"]
                    c = Dict{String,Any}(k => read_attribute(g, k) for k in CASE_ATTRIBUTES)
                    for k in dataset_keys(family)
                        require(haskey(g, k), "$id/$grid/$(k)欠落")
                        c[k] = read(g[k])
                        k != "save_step_index" &&
                            require(read_attribute(g[k], "units") == "1", "units不一致")
                        k in field_keys(family) && require(
                            read_attribute(g[k], "axis_order") == "time,y,x",
                            "場の軸不一致",
                        )
                        if endswith(k, "heat_integral")
                            require(
                                read_attribute(g[k], "axis_order") == "interval,side" &&
                                    read_attribute(g[k], "side_order") ==
                                    "west,east,south,north" &&
                                    read_attribute(g[k], "positive_direction") == "inward",
                                "熱輸送の軸・辺順・符号が不一致",
                            )
                        end
                    end
                    validate_case(c, id, grid, family)
                    grids[grid] = c
                end
                cases[id] = grids
            end
            if official
                expected = expected_cases(family)
                require(
                    Set(keys(cases)) == Set(keys(expected)) &&
                        all(Set(keys(cases[id])) == expected[id] for id in keys(expected)),
                    "公式ケース・格子が不足しています",
                )
            end
            (; metadata = meta, cases)
        end
    catch e
        error("HDF5読取り失敗 $path: $(sprint(showerror,e))")
    end
end

"""
    read_pair(dir)

温度とBurgersの保存場が同じ実行由来か照合する。

# 引数

- `dir`: 保存場を読むディレクトリ。

# 返り値

temperatureとburgersのNamedTuple。run_id・コードhash・Git HEAD混在は拒否する。
"""
function read_pair(dir)
    t = read_fields(joinpath(dir, "temperature.h5"))
    b = read_fields(joinpath(dir, "burgers.h5"))

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        t.metadata["equation_family"] == "temperature" &&
            b.metadata["equation_family"] == "burgers",
        "HDF5の方程式系が不一致",
    )
    for k in ("run_id", "code_sha256", "git_head")
        require(t.metadata[k] == b.metadata[k], "HDF5の計算由来が混在しています: $k")
    end
    (; temperature = t, burgers = b)
end

"""
    input_hashes(dir)

2つの保存場のhashをまとめる。

# 引数

- `dir`: 保存場を読むディレクトリ。

# 返り値

ファイル名からSHA256への辞書。
"""
function input_hashes(dir)
    Dict(n => file_sha(joinpath(dir, n)) for n in ("temperature.h5", "burgers.h5"))
end

"""
    diagnostics(pair, budget)

保存場と辺別積分から熱収支・解析解誤差・収束を評価する。

# 引数

- `pair`: 同じrun_idとコード来歴の温度/Burgers保存場。
- `budget`: 熱量と辺別区間積分から熱収支を返す学生関数。

# 返り値

schema・run_id・来歴・ケース診断・収束を持つ辞書。非有限診断は保存しない。
"""
function diagnostics(pair, budget)
    cases = Dict{String,Any}()
    convergence = Dict{String,Any}()
    for (family, data) in (("temperature", pair.temperature), ("burgers", pair.burgers))
        for (id, grids) in data.cases
            ds = Dict{String,Any}()
            for (grid, c) in grids
                d = Dict{String,Any}(
                    "conditions" => Dict(k => c[k] for k in CASE_ATTRIBUTES),
                    "time" => c["time"],
                    "step_count" => length(c["step_dt"]),
                    "dt_min" => minimum(c["step_dt"]),
                    "dt_max" => maximum(c["step_dt"]),
                )
                if family == "temperature"
                    U = c["temperature"]
                    d["heat"] = [
                        c["dx"] * c["dy"] * sum(@view U[:, :, k]) for
                        k in eachindex(c["time"])
                    ]
                    for (key, fun) in (("minimum", minimum), ("maximum", maximum))
                        d[key] = [fun(@view U[:, :, k]) for k in eachindex(c["time"])]
                    end
                    adv = c["advective_heat_integral"]
                    diff = c["diffusive_heat_integral"]
                    b = budget(c["time"], d["heat"], adv, diff)
                    for (k, v) in pairs(b)
                        require(all(isfinite, v), "$id: 非有限な熱収支診断です")
                        d[string(k)] = v
                    end
                    d["advective_heat_integral"] =
                        Dict(side => collect(adv[j, :]) for (j, side) in enumerate(SIDES))
                    d["diffusive_heat_integral"] =
                        Dict(side => collect(diff[j, :]) for (j, side) in enumerate(SIDES))
                    if id in PERIODIC_CASES ||
                       startswith(id, "comparison_") ||
                       id == "closed_fixed"
                        d["l2_error"] = [
                            sqrt(
                                sum(
                                    abs2,
                                    U[:, :, k] - exact_temperature(c["x"], c["y"], t, id),
                                ) / (c["nx"] * c["ny"]),
                            ) for (k, t) in enumerate(c["time"])
                        ]
                    end
                else
                    for component in ("u", "v")
                        U = c[component]
                        d[component * "_minimum"] =
                            [minimum(@view U[:, :, k]) for k in eachindex(c["time"])]
                        d[component * "_maximum"] =
                            [maximum(@view U[:, :, k]) for k in eachindex(c["time"])]
                        d[component * "_l2_error"] = Float64[]
                    end
                    for (k, t) in enumerate(c["time"])
                        ue, ve = exact_burgers(c["x"], c["y"], t)
                        for (key, exact) in (("u", ue), ("v", ve))
                            push!(
                                d[key * "_l2_error"],
                                sqrt(sum(abs2, c[key][:, :, k] - exact) / length(exact)),
                            )
                        end
                    end
                    d["l2_error"] = hypot.(d["u_l2_error"], d["v_l2_error"])
                end
                ds[grid] = d
            end
            cases[id] = ds
            if id in PERIODIC_CASES || id == "periodic_burgers"
                ids = sort(collect(keys(ds)); by = g->ds[g]["conditions"]["nx"])
                errors = [last(ds[g]["l2_error"]) for g in ids]
                require(all(v->isfinite(v) && v > 0, errors), "最終誤差が不正です")
                convergence[id] = Dict(
                    "grid_ids" => ids,
                    "errors" => errors,
                    "orders" => log2.(errors[1:(end - 1)] ./ errors[2:end]),
                )
            end
        end
    end
    Dict{String,Any}(
        "schema_version" => 1,
        "task_id" => "N07",
        "run_id" => pair.temperature.metadata["run_id"],
        "diagnostics_complete" => true,
        "provenance" => pair.temperature.metadata,
        "cases" => cases,
        "convergence" => convergence,
    )
end

"""
    read_summary(path, input_dir)

診断TOMLのケースと出自が現在の両HDF5と一致するか確認する。

# 引数

- `path`: 対象ファイルのパス。
- `input_dir`: 保存場HDF5を読むディレクトリ。

# 返り値

解析辞書。古い出自や未完成の熱収支は再解析を求めて拒否する。
"""
function read_summary(path, input_dir)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(isfile(path), "summary欠落。analyze.jlを実行してください")
    s = TOML.parsefile(path)

    # 保存済みの値と出自を検証してから、解析や作図に使う。
    data = read_pair(input_dir)

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        get(s, "schema_version", 0) == 1 &&
            get(s, "task_id", "") == "N07" &&
            get(s, "source_sha256", nothing) == input_hashes(input_dir) &&
            get(s, "run_id", nothing) == data.temperature.metadata["run_id"],
        "summaryの出自がHDF5と不一致です。analyze.jlを再実行してください",
    )
    require(get(s, "diagnostics_complete", false)===true, "熱収支解析が未完了です")
    require(
        haskey(s, "cases") &&
            Set(keys(s["cases"])) ==
            union(Set(keys(data.temperature.cases)), Set(keys(data.burgers.cases))),
        "summaryのケースが不足しています",
    )
    s
end

"""
    check_complete(output_dir = DEFAULT_OUTPUT_DIR)

熱収支・解析解誤差・格子収束・図hash・全8出力を確認する。

# 引数

- `output_dir`: 公式成果物を書き出すディレクトリ。

# 返り値

完了時true。小さい熱収支残差だけでは完了とせず、保存場から独立に診断する。
"""
function check_complete(output_dir = DEFAULT_OUTPUT_DIR)
    # 保存済みの値と出自を検証してから、解析や作図に使う。
    pair = read_pair(output_dir)
    s = read_summary(joinpath(output_dir, "summary.toml"), output_dir)
    for (id, conv) in s["convergence"]
        expected = grids_for(id)
        require(length(conv["errors"]) == length(expected) == 3, "3格子が必要です")
        errors = conv["errors"]
        p = conv["orders"]
        lo, hi = id == "periodic_diffusion" ? (1.8, 2.2) : (0.8, 1.2)
        require(
            all(diff(errors) .< 0) &&
                length(p) == 2 &&
                all(lo .<= p .<= hi) &&
                isapprox(p, log2.(errors[1:2] ./ errors[2:3]); rtol = 1e-12),
            "$id: 収束基準を満たしません",
        )
    end

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        Set(keys(s["convergence"])) == Set((PERIODIC_CASES..., "periodic_burgers")),
        "収束条件が不足しています",
    )
    # Verify student diagnostics against persisted inputs, without rerunning a solver.
    for (id, grids) in pair.temperature.cases, (grid, c) in grids
        d = s["cases"][id][grid]
        q = [
            c["dx"] * c["dy"] * sum(@view c["temperature"][:, :, k]) for
            k in eachindex(c["time"])
        ]
        transport =
            sum(abs, c["advective_heat_integral"]) + sum(abs, c["diffusive_heat_integral"])
        tol = 1e-12 * max(1, abs(first(q)), transport)
        require(
            isapprox(d["heat"], q; rtol = 1e-13, atol = 1e-14),
            "$id: 熱量が保存場と不一致です",
        )
        require(
            length(d["residual"]) == length(q) &&
                maximum(abs, d["residual"]) <= tol &&
                isapprox(
                    q .- first(q),
                    d["cumulative_input"] + d["residual"];
                    atol = tol,
                    rtol = 0,
                ),
            "$id: 熱収支が不一致です",
        )
        require(
            first(d["cumulative_input"]) == 0 && isapprox(
                diff(d["cumulative_input"]),
                d["net_input"];
                rtol = 1e-12,
                atol = tol,
            ),
            "$id: 区間と累積収支が不一致です",
        )
        for k in eachindex(d["net_input"])
            require(
                isapprox(
                    d["net_input"][k],
                    sum(c["advective_heat_integral"][:, k]) +
                    sum(c["diffusive_heat_integral"][:, k]);
                    atol = tol,
                    rtol = 0,
                ),
                "$id: 境界輸送と不一致です",
            )
        end
    end
    # Recompute the analytic error from stored fields; a small heat residual alone is insufficient.
    for id in (PERIODIC_CASES..., "periodic_burgers")
        data = id == "periodic_burgers" ? pair.burgers : pair.temperature
        conv = s["convergence"][id]
        for (k, grid) in enumerate(conv["grid_ids"])
            c = data.cases[id][grid]
            if id == "periodic_burgers"
                u, v = exact_burgers(c["x"], c["y"], 1.)
                e = sqrt(
                    (sum(abs2, c["u"][:, :, end] - u) + sum(abs2, c["v"][:, :, end] - v)) / length(u),
                )
            else
                u = exact_temperature(c["x"], c["y"], 1., id)
                e = sqrt(sum(abs2, c["temperature"][:, :, end] - u) / length(u))
            end
            require(
                isapprox(conv["errors"][k], e; rtol = 1e-12),
                "$id: 解析解誤差が保存場と不一致です",
            )
        end
    end
    path = joinpath(output_dir, "plots.toml")
    require(isfile(path), "plots.toml欠落")
    p = TOML.parsefile(path)
    require(
        p["source_sha256"] == input_hashes(output_dir) &&
            p["source_summary_sha256"] == file_sha(joinpath(output_dir, "summary.toml")),
        "図の出自が不一致です。plot.jlを再実行してください",
    )
    require(
        all(n->isfile(joinpath(output_dir, n)), OUTPUT_NAMES),
        "公式8出力が不足しています",
    )
    require(
        all(
            file_sha(joinpath(output_dir, n)) == p["figure_sha256"][n] for
            n in keys(p["figure_sha256"])
        ) &&
            Set(keys(p["figure_sha256"])) ==
            Set(n for n in OUTPUT_NAMES if endswith(n, ".png")),
        "図が記録と不一致です",
    )
    true
end
