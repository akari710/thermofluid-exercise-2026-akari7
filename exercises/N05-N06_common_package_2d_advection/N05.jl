# N05/N06: 抽出前の基準作成 → 独立プロセスでの再計算 → 数値回帰の比較。基準は上書きしません。

module N05Regression
using TOML, SHA
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "results", "N05")
const DIRECTORIES = (
    "N01_linear_advection",
    "N02_nonlinear_advection",
    "N03_diffusion",
    "N04_advection_diffusion",
)

"""
    sha(path)

ファイル内容のSHA256を求める。

# 引数

- `path`: 対象ファイルのパス。

# 返り値

16進数のString。
"""
sha(path) = bytes2hex(sha256(read(path)))

"""
    plain(x::NamedTuple)

計算結果をTOMLで保存できる値へ再帰的に変換する。

# 引数

- `x`: TOMLへ保存する計算結果。入力は変更しない。

# 返り値

NamedTupleは文字列キーの辞書、配列は要素ごとの変換、SymbolはString。その他の値はそのまま返す。
"""
plain(x::NamedTuple) = Dict(string(k) => plain(v) for (k, v) in pairs(x))

"""
    plain(x::AbstractArray)

計算結果をTOMLで保存できる値へ再帰的に変換する。

# 引数

- `x`: TOMLへ保存する計算結果。入力は変更しない。

# 返り値

NamedTupleは文字列キーの辞書、配列は要素ごとの変換、SymbolはString。その他の値はそのまま返す。
"""
plain(x::AbstractArray) = map(plain, x)

"""
    plain(x::Symbol)

計算結果をTOMLで保存できる値へ再帰的に変換する。

# 引数

- `x`: TOMLへ保存する計算結果。入力は変更しない。

# 返り値

NamedTupleは文字列キーの辞書、配列は要素ごとの変換、SymbolはString。その他の値はそのまま返す。
"""
plain(x::Symbol) = string(x)

"""
    plain(x)

計算結果をTOMLで保存できる値へ再帰的に変換する。

# 引数

- `x`: TOMLへ保存する計算結果。入力は変更しない。

# 返り値

NamedTupleは文字列キーの辞書、配列は要素ごとの変換、SymbolはString。その他の値はそのまま返す。
"""
plain(x) = x

"""
    code_hashes(root)

srcとN01〜N04のJuliaソースの来歴を収集する。

# 引数

- `root`: 学生リポジトリのルート。

# 返り値

リポジトリ相対パスからSHA256への辞書。
"""
function code_hashes(root)
    paths = String[]
    for dir in
        vcat([joinpath(root, "src")], [joinpath(root, "exercises", d) for d in DIRECTORIES])
        for (base, _, names) in walkdir(dir), name in names
            endswith(name, ".jl") && push!(paths, joinpath(base, name))
        end
    end
    Dict(relpath(p, root) => sha(p) for p in sort(paths))
end

"""
    capture(root)

N01〜N04を独立したモジュールへ読み込み、基準計算を行う。

# 引数

- `root`: 学生リポジトリのルート。

# 返り値

schema_version・選択モデル・Git HEAD・コードhash・計算結果の辞書。
"""
function capture(root)
    sandbox = Module(gensym(:Baseline))
    for dir in DIRECTORIES
        Base.include(sandbox, joinpath(root, "exercises", dir, "run.jl"))
    end
    Base.invokelatest(capture_loaded, root, sandbox)
end

"""
    capture_loaded(root, sandbox)

読み込み済みの4課題で公式条件と追加条件の数値結果を集める。

# 引数

- `root`: 学生リポジトリのルート。
- `sandbox`: 対象課題を読み込んだ独立したモジュール。

# 返り値

基準または回帰比較に使う結果と来歴の辞書。
"""
function capture_loaded(root, sandbox)
    a, b, c, d=(
        getfield(sandbox, n) for n in
        (:N01LinearAdvection, :N02NonlinearAdvection, :N03Diffusion, :N04AdvectionDiffusion)
    )
    # Called in a new process; invokelatest bridges included module definitions.
    call(f; kwargs...) = Base.invokelatest(f; kwargs...)
    model = d.SELECTED_MODEL
    results = Dict{String,Any}()
    for scheme in (:upwind, :centered)
        results["N01_$scheme"] = plain(call(a.simulate; scheme))
    end
    for boundary in (:fixed, :periodic)
        results["N02_$boundary"] =
            plain(call(b.simulate; boundary, nx = boundary == :fixed ? 81 : 80))
    end
    results["N02_high_cfl"] =
        plain(call(b.simulate; boundary = :periodic, nx = 40, cfl = 1.1, t_final = 0.02))
    for boundary in (:fixed, :insulated), initial in (:pulse, :mode)
        results["N03_$(boundary)_$initial"] = plain(call(c.simulate; boundary, initial))
    end
    results["N03_high_fo"] = plain(call(c.simulate; fo = 0.6, t_final = 0.01))
    combined = call(d.simulate; model)
    results["N04_combined"] = plain(combined)
    results["N04_advection"] =
        plain(call(d.simulate; model, diffusivity = 0., dt = combined.dt))
    results["N04_diffusion"] =
        plain(call(d.simulate; model, advection = false, dt = combined.dt))
    results["N04_mode"] = plain(call(d.simulate; model, initial = :mode))
    head = try
        strip(read(`git -C $root rev-parse HEAD`, String))
    catch
        "unavailable"
    end
    Dict(
        "schema_version" => 1,
        "selected_model" => string(model),
        "git_head" => head,
        "code_sha256" => code_hashes(root),
        "results" => results,
    )
end

"""
    fresh_capture(root)

新しいJuliaプロセスで数値結果を再計算する。

# 引数

- `root`: 学生リポジトリのルート。

# 返り値

一時TOMLから読んだ結果辞書。子プロセスの失敗は伝わり、公式成果物は変更しない。
"""
function fresh_capture(root)
    mktempdir() do tmp
        path = joinpath(tmp, "capture.toml")
        command = `$(Base.julia_cmd()) --startup-file=no --project=$root $(@__FILE__) _capture $root $path`
        run(command)
        TOML.parsefile(path)
    end
end

"""
    compare!(differences, a, b, path = "results")

辞書・配列を再帰的に比較し、浮動小数の絶対差を記録する。

# 引数

- `differences`: 絶対差を追記する辞書。
- `a`: 保存済みの基準値。
- `b`: 今回計算した比較対象の値。
- `path`: 差分を記録する項目名の階層。ファイルパスではない。

# 返り値

辞書・配列ではnothing、浮動小数では記録した絶対差、その他の一致ではtrue。differencesへ差を追記する。キー・形状・型・数値の不一致はエラー。
"""
function compare!(differences, a, b, path = "results")
    if a isa AbstractDict
        b isa AbstractDict && Set(keys(a)) == Set(keys(b)) || error("N05 キー不一致: $path")
        for k in sort(collect(keys(a)))
            compare!(differences, a[k], b[k], "$path.$k")
        end
    elseif a isa AbstractArray
        b isa AbstractArray && size(a) == size(b) || error("N05 形状不一致: $path")
        for i in eachindex(a)
            compare!(differences, a[i], b[i], "$path[$i]")
        end
    elseif a isa AbstractFloat
        b isa AbstractFloat &&
        isfinite(a) &&
        isfinite(b) &&
        isapprox(a, b; rtol = 1e-12, atol = 1e-13) ||
            error("N05 数値不一致: $path ($a != $b)")
        differences[path] = abs(a - b)
    else
        typeof(a) == typeof(b) && a == b || error("N05 値不一致: $path ($a != $b)")
    end
end

"""
    write_toml(path, data)

同じ出力先の一時ファイルを完成させてTOMLを置き換える。

# 引数

- `path`: 対象ファイルのパス。
- `data`: 検証済みの保存場とmetadata。

# 返り値

nothing。pathを書き換える。生成失敗時は公式パスへの移動を行わず、一時ファイルを片付ける。
"""
function write_toml(path, data)
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        TOML.print(io, data; sorted = true)
        close(io)
        mv(temporary, path; force = true)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

"""
    baseline(; root = ROOT, output_dir = DEFAULT_OUTPUT_DIR)

抽出前の結果とコード来歴を、上書き禁止の基準として保存する。

# 引数

- `root`: 学生リポジトリのルート。
- `output_dir`: 公式成果物を書き出すディレクトリ。

# 返り値

baseline.tomlのパス。既存基準があればエラーにする。
"""
function baseline(; root = ROOT, output_dir = DEFAULT_OUTPUT_DIR)
    path = joinpath(output_dir, "baseline.toml")
    ispath(path) && error("N05 baselineは上書きできません: $path")
    data = fresh_capture(root)
    ispath(path) && error("N05 baselineが既に存在します: $path")
    write_toml(path, data)
    path
end

"""
    verify(;
        root = ROOT,
        baseline_path = joinpath(DEFAULT_OUTPUT_DIR, "baseline.toml"),
        output_dir = DEFAULT_OUTPUT_DIR,
    )

現在の数値結果を抽出前の基準と比較する。

# 引数

- `root`: 学生リポジトリのルート。
- `baseline_path`: 抽出前の数値結果を保存したbaseline TOMLのパス。
- `output_dir`: 公式成果物を書き出すディレクトリ。

# 返り値

回帰結果の辞書。regression.tomlを保存する。基準欠落・破損・schemaや選択モデルの変更・数値不一致では保存しない。
"""
function verify(;
    root = ROOT,
    baseline_path = joinpath(DEFAULT_OUTPUT_DIR, "baseline.toml"),
    output_dir = DEFAULT_OUTPUT_DIR,
)
    isfile(baseline_path) || error(
        "N05 baseline欠落: $(baseline_path)。抽出前commitの隔離コピーでN05.jl baselineを実行してください。",
    )
    before = try
        TOML.parsefile(baseline_path)
    catch e
        error("N05 baseline破損: $baseline_path: $e")
    end
    get(before, "schema_version", 0) == 1 ||
        error("N05 baseline schema不一致: $baseline_path")
    now = fresh_capture(root)
    get(before, "selected_model", nothing) == now["selected_model"] ||
        error("N05 N04選択変更: $baseline_path")
    differences = Dict{String,Float64}()
    compare!(differences, before["results"], now["results"])

    # 計算値と条件を、保存する診断情報にまとめる。
    report = Dict(
        "schema_version" => 1,
        "baseline_sha256" => sha(baseline_path),
        "code_sha256" => now["code_sha256"],
        "git_head" => now["git_head"],
        "selected_model" => now["selected_model"],
        "differences" => differences,
    )
    write_toml(joinpath(output_dir, "regression.toml"), report)
    report
end

"""
    main(args = ARGS)

引数に応じて基準作成、回帰確認、内部captureを選ぶ。

# 引数

- `args`: 入口へ渡す引数の列。

# 返り値

選択した処理の結果。引数なしはverify。不明な引数は使い方を示して停止する。
"""
function main(args = ARGS)
    isempty(args) && return verify()
    args == ["baseline"] && return baseline()
    args == ["verify"] && return verify()
    if length(args) == 3 && args[1] == "_capture"
        return write_toml(args[3], capture(args[2]))
    end
    error("使い方: N05.jl [baseline|verify]")
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
