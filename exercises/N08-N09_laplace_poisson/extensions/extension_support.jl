# N08/N09: 発展問題と独立診断を提供し、必修とは別のschema・出力を検査します。

include("../provided_support.jl")
const NEUMANN_CASES = (
    "mixed_laplace_flux",
    "mixed_poisson_flux",
    "mixed_poisson_zero_flux",
    "pure_poisson_flux",
    "pure_poisson_zero_flux",
)

"""
    trapezoid_weights(nx, ny, dx, dy)

半セル境界を持つ台形積分の節点重みを作る。

# 引数

- `nx`: x方向の格子点数。
- `ny`: y方向の格子点数。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。

# 返り値

(nx,ny)の重み行列。
"""
function trapezoid_weights(nx, ny, dx, dy)
    a = ones(nx)
    b = ones(ny)
    a[[1, end]].=0.5
    b[[1, end]].=0.5
    dx * dy .* (a * b')
end

"""
    weighted_mean(u, dx, dy)

節点の半セル重みを使って場の平均を求める。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。

# 返り値

重み付き平均のスカラー値。uは保持する。
"""
weighted_mean(u, dx, dy) =
    sum(trapezoid_weights(size(u)..., dx, dy) .* u) /
    sum(trapezoid_weights(size(u)..., dx, dy))

"""
    neumann_problem(id, nx, ny)

混合または純Neumannの独立した問題と解析解を作る。

# 引数

- `id`: 格子点数と対応するケース名。
- `nx`: x方向の格子点数。
- `ny`: y方向の格子点数。

# 返り値

x,y,dx,dy,exact,f,bc,u0のNamedTuple。Neumann値は外向き微分。
"""
function neumann_problem(id, nx, ny)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(id in NEUMANN_CASES, "case_id不正")

    # 座標と格子幅を用意し、配列の添字との対応をそろえる。
    x = collect(range(0., 2.; length = nx))
    y = collect(range(0., 1.; length = ny))
    dx = 2 / (nx - 1)
    dy = 1 / (ny - 1)
    k = pi / 2
    d(v) = (; kind = :dirichlet, values = Float64.(v))
    n(v) = (; kind = :neumann, values = Float64.(v))
    if startswith(id, "mixed")
        poisson = id == "mixed_poisson_flux"
        exact = [H(a, b) + (poisson ? V(a, b) : 0.) for a in x, b in y]
        f =
            poisson ? [-((pi / 2) ^ 2 + pi ^ 2) * V(a, b) for a in x, b in y] :
            zeros(nx, ny)
        south =
            [-1 - sin(k * a) * k / sinh(k) - (poisson ? pi * sin(k * a) : 0.) for a in x]
        north = [
            1 + sin(k * a) * k * cosh(k) / sinh(k) - (poisson ? pi * sin(k * a) : 0.)
            for a in x
        ]
        if id == "mixed_poisson_zero_flux"
            exact = [a ^ 2 for a in x, b in y]
            f = fill(2., nx, ny)
            south = zeros(nx)
            north = zeros(nx)
        end
        bc = (;
            west = d(exact[1, :]),
            east = d(exact[end, :]),
            south = n(south),
            north = n(north),
        )
    elseif id == "pure_poisson_flux"
        exact = [a ^ 2 + 2b ^ 2 for a in x, b in y]
        f = fill(6., nx, ny)
        bc = (;
            west = n(zeros(ny)),
            east = n(fill(4., ny)),
            south = n(zeros(nx)),
            north = n(fill(4., nx)),
        )
    else
        exact = [cos(pi * a) + 0.5cos(2pi * b) for a in x, b in y]
        f = [-pi ^ 2 * cos(pi * a) - 2pi ^ 2 * cos(2pi * b) for a in x, b in y]
        bc = (;
            west = n(zeros(ny)),
            east = n(zeros(ny)),
            south = n(zeros(nx)),
            north = n(zeros(nx)),
        )
    end
    (; x, y, dx, dy, exact, f, bc, u0 = zeros(nx, ny))
end

"""
    is_pure(bc)

全ての辺がNeumann条件か判定する。

# 引数

- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとvaluesを持ち、Neumann値は外向き法線微分。

# 返り値

Bool。
"""
is_pure(bc) = all(b->b.kind == :neumann, values(bc))

"""
    independent_compatibility(f, dx, dy, bc)

右辺の領域積分と外向き境界流束の不一致を独立に求める。

# 引数

- `f`: Poisson方程式の右辺行列。入力は変更しない。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとvaluesを持ち、Neumann値は外向き法線微分。

# 返り値

defectと判定用scaleのNamedTuple。
"""
function independent_compatibility(f, dx, dy, bc)
    nx, ny = size(f)
    wx = ones(nx)
    wy = ones(ny)
    wx[[1, end]].=0.5
    wy[[1, end]].=0.5
    edges = (
        dy .* wy .* bc.west.values,
        dy .* wy .* bc.east.values,
        dx .* wx .* bc.south.values,
        dx .* wx .* bc.north.values,
    )
    defect = sum(trapezoid_weights(nx, ny, dx, dy) .* f) - sum(sum(e) for e in edges)
    (;
        defect,
        scale = max(
            1.,
            sum(abs, trapezoid_weights(nx, ny, dx, dy) .* f),
            sum(sum(abs, e) for e in edges),
        ),
    )
end

"""
    extension_residual(u, f, dx, dy, bc)

発展の保存場から境界を含む残差を独立に再計算する。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `f`: Poisson方程式の右辺行列。入力は変更しない。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。
- `bc`: west,east,south,north順の境界NamedTuple。各辺はkindとvaluesを持ち、Neumann値は外向き法線微分。

# 返り値

固定点を除く最大絶対残差。
"""
function extension_residual(u, f, dx, dy, bc)
    nx, ny = size(u)
    norm = 0.
    for j in 1:ny, i in 1:nx
        fixed =
            (i == 1 && bc.west.kind == :dirichlet)||(i == nx && bc.east.kind == :dirichlet)||(
                j == 1 && bc.south.kind == :dirichlet
            )||(j == ny && bc.north.kind == :dirichlet)
        fixed && continue
        l = i == 1 ? u[2, j] + 2dx * bc.west.values[j] : u[i - 1, j]
        r = i == nx ? u[nx - 1, j] + 2dx * bc.east.values[j] : u[i + 1, j]
        b = j == 1 ? u[i, 2] + 2dy * bc.south.values[i] : u[i, j - 1]
        t = j == ny ? u[i, ny - 1] + 2dy * bc.north.values[i] : u[i, j + 1]
        norm = max(
            norm,
            abs(f[i, j] - (l - 2u[i, j] + r) / dx ^ 2 - (b - 2u[i, j] + t) / dy ^ 2),
        )
    end
    norm
end

"""
    normal_derivatives(u, dx, dy)

各辺の外向き法線微分を片側差分で評価する。

# 引数

- `u`: 場の値。配列の添字は座標の順に対応する。
- `dx`: x方向の有限な正の格子幅。
- `dy`: y方向の有限な正の格子幅。

# 返り値

west,east,south,northのベクトルを持つNamedTuple。
"""
normal_derivatives(u, dx, dy) = (;
    west = (-3u[1, :] + 4u[2, :] - u[3, :]) ./ (-2dx),
    east = (3u[end, :] - 4u[end - 1, :] + u[end - 2, :]) ./ (2dx),
    south = (-3u[:, 1] + 4u[:, 2] - u[:, 3]) ./ (-2dy),
    north = (3u[:, end] - 4u[:, end - 1] + u[:, end - 2]) ./ (2dy),
)

"""
    extension_case(
        p,
        result,
        id;
        method = "weighted_jacobi",
        omega = 2 / 3,
        elapsed = 0.,
    )

収束済みの発展場と方法・時間・履歴を保存ケースにまとめる。

# 引数

- `p`: 座標・格子幅・解析解・右辺と境界を持つ問題データ。
- `result`: シミュレーションの結果と診断量。
- `id`: 格子点数と対応するケース名。
- `method`: 反復法の識別子。
- `omega`: 更新の緩和係数。許容範囲は下記の検証に従う。
- `elapsed`: warm-upを除く計算時間（秒）。

# 返り値

ケース辞書。未収束は拒否する。Neumannでは境界・適合条件と重み付き平均も記録する。
"""
function extension_case(
    p,
    result,
    id;
    method = "weighted_jacobi",
    omega = 2 / 3,
    elapsed = 0.,
)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(result.converged, "発展未収束: $id $(size(p.f)) $(result.reason)")
    c = Dict{String,Any}(
        "case_id" => id,
        "nx" => length(p.x),
        "ny" => length(p.y),
        "dx" => p.dx,
        "dy" => p.dy,
        "converged" => result.converged,
        "reason" => string(result.reason),
        "iterations" => result.iterations,
        "threshold" => result.threshold,
        "initial_residual" => first(result.residual_history),
        "x" => p.x,
        "y" => p.y,
        "u" => result.u,
        "f" => p.f,
        "g" => hasproperty(p, :g) ? p.g : zeros(size(p.f)),
        "iteration" => collect(0:result.iterations),
        "residual_history" => result.residual_history,
        "update_history" => result.update_history,
        "method" => method,
        "omega" => omega,
        "elapsed_seconds" => elapsed,
    )
    if hasproperty(p, :bc)
        c["bc"] = p.bc
        c["compatibility_defect"] =
            is_pure(p.bc) ? independent_compatibility(p.f, p.dx, p.dy, p.bc).defect : 0.
        c["weighted_mean"] = weighted_mean(result.u, p.dx, p.dy)
    end
    c
end

"""
    extension_metadata(ext, run_id)

必修と分離した発展ID・手法とコード来歴をまとめる。

# 引数

- `ext`: 発展課題の識別子（neumannまたはrelaxation）。
- `run_id`: 同じ実行の保存場を対応付ける識別子。

# 返り値

ルート属性の辞書。
"""
function extension_metadata(ext, run_id)
    m = source_metadata("N09", run_id)
    m["extension_id"] = ext
    m["boundary"] = ext == "neumann" ? "per_case" : "dirichlet"
    m["method"] = ext == "neumann" ? "weighted_jacobi" : "per_case"
    hashes = TOML.parse(m["code_sha256"])
    for name in ("extension_support.jl", ext * ".jl", "run_" * ext * ".jl")
        hashes["extensions/" * name] = file_sha(joinpath(@__DIR__, name))
    end
    m["code_sha256"] = sprint(io->TOML.print(io, hashes; sorted = true))
    m
end

"""
    validate_extension_case(c, meta)

発展の場・境界・適合条件・平均0・残差と履歴を検証する。

# 引数

- `c`: 格子・座標・条件・保存場を持つケース辞書。
- `meta`: HDF5のルート属性。

# 返り値

L2/Linf誤差・独立残差・反復数・時間・停止理由の辞書。不整合や未収束は拒否する。
"""
function validate_extension_case(c, meta)
    ext = meta["extension_id"]
    nx, ny = c["nx"], c["ny"]

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require((nx, ny) in OFFICIAL_GRIDS, "発展格子不一致")
    p = ext == "neumann" ? neumann_problem(c["case_id"], nx, ny) : problem("N09", nx, ny)
    require(
        c["x"] == p.x && c["y"] == p.y && c["dx"] == p.dx && c["dy"] == p.dy,
        "発展座標不一致",
    )
    for k in ("u", "f", "g")
        require(
            c[k] isa Matrix{Float64} && size(c[k]) == (nx, ny) && all(isfinite, c[k]),
            "発展場不正: $k",
        )
    end
    require(
        c["f"] == p.f && c["g"] == (ext == "neumann" ? zeros(nx, ny) : p.g),
        "発展右辺／境界不一致",
    )
    n = c["iterations"]
    rh = c["residual_history"]
    uh = c["update_history"]
    require(
        n isa Integer &&
            0 <= n <= meta["maxiter"] &&
            c["iteration"] == collect(0:n) &&
            length(rh) == length(uh) == n + 1,
        "発展履歴長不一致",
    )
    require(
        rh isa Vector{Float64} &&
            uh isa Vector{Float64} &&
            all(v->isfinite(v) && v >= 0, rh) &&
            all(v->isfinite(v) && v >= 0, uh) &&
            first(uh) == 0.,
        "発展履歴不正",
    )
    require(
        c["converged"] == true &&
            c["reason"] == "converged" &&
            c["initial_residual"] == first(rh) &&
            c["threshold"] == max(meta["atol"], meta["rtol"] * first(rh)) &&
            last(rh) <= c["threshold"],
        "発展停止不正",
    )
    require(isfinite(c["elapsed_seconds"]) && c["elapsed_seconds"] >= 0, "計測時間不正")
    exact = copy(p.exact)
    if ext == "neumann"
        require(
            c["bc"] == p.bc && c["method"] == "weighted_jacobi" && c["omega"] == 2 / 3,
            "発展境界／omega不一致",
        )
        if is_pure(p.bc)
            compat = independent_compatibility(c["f"], p.dx, p.dy, p.bc)
            require(
                abs(compat.defect) <= 1e-12 * compat.scale &&
                    c["compatibility_defect"] == compat.defect,
                "保存可解条件違反",
            )
            require(abs(weighted_mean(c["u"], p.dx, p.dy)) < 1e-13, "保存平均0不一致")
            exact.-=weighted_mean(exact, p.dx, p.dy)
        end
        require(
            c["weighted_mean"] == weighted_mean(c["u"], p.dx, p.dy),
            "保存平均情報不一致",
        )
        residual = extension_residual(c["u"], p.f, p.dx, p.dy, p.bc)
        for side in (:west, :east)
            values = side == :west ? c["u"][1, :] : c["u"][end, :]
            p.bc[side].kind == :dirichlet &&
                require(values == p.bc[side].values, "混合Dirichlet辺不一致")
        end
        err = c["u"] - exact
    else
        require(
            c["case_id"] in ("jacobi", "gauss_seidel", "sor") &&
                c["method"] == c["case_id"],
            "緩和method不一致",
        )
        require(c["omega"] == (c["case_id"] == "sor" ? 1.5 : 1.), "緩和omega不一致")
        require(
            c["u"][1, :] == p.g[1, :] &&
                c["u"][end, :] == p.g[end, :] &&
                c["u"][:, 1] == p.g[:, 1] &&
                c["u"][:, end] == p.g[:, end],
            "発展Dirichlet辺不一致",
        )
        residual = independent_residual(c["u"], p.f, p.dx, p.dy)
        err = c["u"][2:(end - 1), 2:(end - 1)] - exact[2:(end - 1), 2:(end - 1)]
    end
    require(
        isapprox(residual, last(rh); atol = 3e-11, rtol = 1e-5) &&
            residual <= c["threshold"] + 3e-11,
        "発展独立残差不一致",
    )
    Dict(
        "l2_error" => sqrt(sum(abs2, err) / length(err)),
        "linf_error" => maximum(abs, err),
        "independent_residual" => residual,
        "iterations" => n,
        "elapsed_seconds" => c["elapsed_seconds"],
        "reason" => c["reason"],
    )
end

"""
    validate_extension_metadata(meta)

発展IDと境界・手法・共通schemaを照合する。

# 引数

- `meta`: HDF5のルート属性。

# 返り値

成功時nothing。保存済みhashは書式を検査し、現在ソースとのbyte一致は要求しない。
"""
function validate_extension_metadata(meta)
    ext = meta["extension_id"]

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(ext in ("neumann", "relaxation"), "発展ID不正")
    expected = extension_metadata(ext, meta["run_id"])
    # Saved provenance remains valid after source edits; its digest syntax is checked separately.
    normal = copy(meta)
    normal["boundary"] = "dirichlet"
    normal["method"] = "jacobi"
    validate_metadata(normal)
    require(
        meta["boundary"] == expected["boundary"] && meta["method"] == expected["method"],
        "発展root不一致",
    )
end

"""
    write_extension_fields(path, cases, meta)

発展ケースを圧縮HDF5に保存して再検証する。

# 引数

- `path`: 対象ファイルのパス。
- `cases`: ケース名をキーにした計算結果の辞書。
- `meta`: HDF5のルート属性。

# 返り値

metadataとcasesのNamedTuple。Julia場は(x,y)、属性axis_orderはy,x。
"""
function write_extension_fields(path, cases, meta)
    # 計算や書込みの前に、入力条件をまとめて確認する。
    validate_extension_metadata(meta)
    for c in values(cases)
        validate_extension_case(c, meta)
    end
    h5open(path, "w") do h
        for (k, v) in meta
            attributes(h)[k] = v
        end
        groups = create_group(h, "cases")
        for (key, c) in cases
            g = create_group(groups, key)
            for (k, v) in c
                if k in DATA_KEYS
                    if v isa AbstractVector
                        g[
                            k,
                            chunk = (min(length(v), 16384),),
                            shuffle = true,
                            deflate = 9,
                        ] = v
                    else
                        g[k, chunk = size(v), shuffle = true, deflate = 9] = v
                    end
                    attributes(g[k])["units"] = k == "iteration" ? "iteration" : "1"
                    k in ("u", "f", "g") && (attributes(g[k])["axis_order"] = "y,x")
                elseif k == "bc"
                    b = create_group(g, "bc")
                    for side in keys(v)
                        edge = create_group(b, string(side))
                        attributes(edge)["kind"] = string(v[side].kind)
                        edge["values"] = v[side].values
                        attributes(edge["values"])["units"] = "1"
                    end
                else
                    attributes(g)[k] = v
                end
            end
        end
    end
    read_extension_fields(path)
end

"""
    read_extension_fields(path)

発展HDF5の軸・単位・全ケース集合と独立診断を検証して読む。

# 引数

- `path`: 対象ファイルのパス。

# 返り値

metadataとcasesのNamedTuple。不正データは拒否する。
"""
function read_extension_fields(path)
    h5open(path, "r") do h
        meta = Dict{String,Any}(
            k => read_attribute(h, k) for k in (ROOT_KEYS..., "extension_id")
        )
        validate_extension_metadata(meta)
        cases = Dict{String,Any}()
        for key in keys(h["cases"])
            g = h["cases/$key"]
            c = Dict{String,Any}(k => read_attribute(g, k) for k in keys(attributes(g)))
            for k in DATA_KEYS
                c[k] = read(g[k])
                require(
                    read_attribute(g[k], "units") == (k == "iteration" ? "iteration" : "1"),
                    "発展units不一致",
                )
                k in ("u", "f", "g") &&
                    require(read_attribute(g[k], "axis_order") == "y,x", "発展軸不一致")
            end
            if meta["extension_id"] == "neumann"
                c["bc"] = NamedTuple{(:west, :east, :south, :north)}(
                    Tuple(
                        (;
                            kind = Symbol(read_attribute(g["bc/" * string(side)], "kind")),
                            values = read(g["bc/" * string(side) * "/values"]),
                        ) for side in (:west, :east, :south, :north)
                    ),
                )
                for side in keys(c["bc"])
                    require(
                        read_attribute(g["bc/" * string(side) * "/values"], "units") == "1",
                        "法線微分単位不一致",
                    )
                end
            end
            validate_extension_case(c, meta)
            cases[key] = c
        end
        ids =
            meta["extension_id"] == "neumann" ? NEUMANN_CASES :
            ("jacobi", "gauss_seidel", "sor")
        require(
            Set(keys(cases)) ==
            Set(id * "_" * case_id(n...) for id in ids for n in OFFICIAL_GRIDS),
            "発展ケース集合不一致",
        )
        (; metadata = meta, cases)
    end
end

"""
    extension_diagnostics(pair, hash)

発展の独立診断と3格子の約2次収束をまとめる。

# 引数

- `pair`: 同じrun_idとコード来歴を持つLaplace/Poissonの保存場。発展課題では選択した発展の保存場。
- `hash`: 解析元fields.h5のSHA256。

# 返り値

metadata・保存場hash・ケース別診断・次数を持つ辞書。
"""
function extension_diagnostics(pair, hash)
    meta = pair.metadata
    d = Dict(k => validate_extension_case(c, meta) for (k, c) in pair.cases)
    orders = Dict{String,Any}()
    ids =
        meta["extension_id"] == "neumann" ? ("pure_poisson_zero_flux",) :
        ("jacobi", "gauss_seidel", "sor")
    for id in ids, norm in ("l2", "linf")
        errs = [d[id * "_" * case_id(n...)][norm * "_error"] for n in OFFICIAL_GRIDS]
        p = log2.(errs[1:2] ./ errs[2:3])
        require(all(v->1.8 <= v <= 2.2, p), "発展3格子次数不一致")
        orders[id * "_" * norm] = p
    end
    merge(
        copy(meta),
        Dict(
            "source_fields_sha256" => hash,
            "diagnostics_complete" => true,
            "cases" => d,
            "orders" => orders,
        ),
    )
end

"""
    check_extension(dir)

発展診断と図の来歴・hashを保存場と照合する。

# 引数

- `dir`: 保存場を読むディレクトリ。

# 返り値

完了時true。
"""
function check_extension(dir)
    pair = read_extension_fields(joinpath(dir, "fields.h5"))
    summary = TOML.parsefile(joinpath(dir, "summary.toml"))
    plots = TOML.parsefile(joinpath(dir, "plots.toml"))

    # 計算や書込みの前に、入力条件をまとめて確認する。
    require(
        summary == extension_diagnostics(pair, file_sha(joinpath(dir, "fields.h5"))),
        "発展summary不一致",
    )
    require(
        plots["run_id"] == summary["run_id"] &&
            plots["source_fields_sha256"] == file_sha(joinpath(dir, "fields.h5")) &&
            plots["source_summary_sha256"] == file_sha(joinpath(dir, "summary.toml")),
        "発展来歴不一致",
    )
    require(
        all(plots["png_sha256"][n] == file_sha(joinpath(dir, n)) for n in PNG_FILES),
        "発展PNG不一致",
    )
    true
end
