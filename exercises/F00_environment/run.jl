# F00の提供コード。環境の観測 → 結果の表示 → 手動確認を含む進捗判定の順に読みます。
# このファイルに数値TODOはありません。直接実行は診断だけを行います。

if !isdefined(Main, :CourseWorkflow)
    include(joinpath(@__DIR__, "..", "..", "scripts", "lib", "CourseWorkflow.jl"))
end

module F00Environment

using Main.CourseWorkflow

export ObservedCheck,
    PreflightReport,
    REQUIRED_JULIA_VERSION,
    SUPPORTED_AGENTS,
    classify_runtime,
    collect_preflight,
    parse_preflight_arguments,
    print_preflight,
    run_f00_preflight

const REQUIRED_JULIA_VERSION = v"1.13.0"
const SUPPORTED_AGENTS = ("copilot", "codex", "amazon-q")

"""
    ObservedCheck(id, passed, observed, action)

1つの環境確認の判定と、受講生に表示する観測値・対応方法を保持する。

# 引数

- `id`: 確認項目の `Symbol`。
- `passed`: 確認に合格したかを表す `Bool`。
- `observed`: 観測値の `String`。
- `action`: 要対応時に表示する `String`。

# 返り値

4つの引数を同名のフィールドに持つ `ObservedCheck`。
"""
struct ObservedCheck
    id::Symbol
    passed::Bool
    observed::String
    action::String
end

"""
    PreflightReport(runtime, workspace, julia, git)

実行環境・作業場所・Julia・Gitの確認結果をまとめる。

# 引数

- `runtime`: 実行環境の `ObservedCheck`。
- `workspace`: 作業場所の `ObservedCheck`。
- `julia`: Juliaの版の `ObservedCheck`。
- `git`: Gitコマンドの `ObservedCheck`。

# 返り値

4つの確認結果を同名のフィールドに持つ `PreflightReport`。
"""
struct PreflightReport
    runtime::ObservedCheck
    workspace::ObservedCheck
    julia::ObservedCheck
    git::ObservedCheck
end

"""
    classify_runtime(platform, kernel_release)

観測したOSとカーネル名から、教材を実行できる環境か分類する。

# 引数

- `platform`: `:windows`、`:macos`、`:linux` などのOS識別子。
- `kernel_release`: カーネル名の文字列。空文字の場合はLinuxを判定不能とする。

# 返り値

`kind`、`passed`、`observed`、`action` を持つ `NamedTuple`。判定できない環境は `passed=false` と対応方法を返す。
"""
function classify_runtime(platform, kernel_release)
    # OSの識別子とカーネル名を正規化してから環境を分類する。
    platform_name = Symbol(platform)
    release = strip(String(kernel_release))
    lower_release = lowercase(release)

    # 同じLinuxでもWSL1・WSL2・native Linuxを区別する。
    if platform_name == :windows
        return (
            kind = :native_windows,
            passed = false,
            observed = "native Windows",
            action = "WindowsではWSL2 Ubuntu 24.04を起動し、Linux側のJulia・Git・SSH・agentを使用してください。",
        )
    elseif platform_name == :macos
        return (kind = :macos, passed = true, observed = "native macOS", action = "")
    elseif platform_name == :linux
        isempty(release) && return (
            kind = :unknown_linux,
            passed = false,
            observed = "Linuxのkernel releaseを判定できません",
            action = "WSL2またはnative Linuxの端末でkernel情報を確認し、F00を再実行してください。",
        )
        if occursin("microsoft", lower_release) && occursin("wsl2", lower_release)
            return (
                kind = :wsl2,
                passed = true,
                observed = "WSL2 Ubuntu 24.04相当 (kernel: $release)",
                action = "",
            )
        elseif occursin("microsoft", lower_release)
            return (
                kind = :wsl1,
                passed = false,
                observed = "WSL1相当 (kernel: $release)",
                action = "WSL2 Ubuntu 24.04へ更新してから、Linux側でF00を再実行してください。",
            )
        else
            return (
                kind = :native_linux,
                passed = true,
                observed = "native Linux (kernel: $release)",
                action = "",
            )
        end
    end

    (
        kind = :unknown,
        passed = false,
        observed = "未対応の実行環境: $platform_name",
        action = "WSL2 Ubuntu 24.04、native macOS、またはnative LinuxでF00を再実行してください。",
    )
end

"""
    default_platform_probe()

実行中のOSをJuliaのシステム情報から取得する。

# 引数

なし。

# 返り値

`:windows`、`:macos`、`:linux`、未対応環境では `:unknown`。
"""
default_platform_probe() =
    Sys.iswindows() ? :windows : Sys.isapple() ? :macos : Sys.islinux() ? :linux : :unknown

"""
    default_kernel_probe()

Linuxのカーネル名を、読取り可能なシステムファイルから取得する。

# 引数

なし。

# 返り値

カーネル名の `String`。Linux以外やファイルがない場合は空文字列。
"""
function default_kernel_probe()
    if Sys.islinux()
        path = "/proc/sys/kernel/osrelease"
        return isfile(path) ? strip(read(path, String)) : ""
    end
    ""
end

"""
    default_runtime_probe()

環境分類と作業場所の確認に使う観測値をまとめる。

# 引数

なし。

# 返り値

`platform`、`kernel_release`、`workspace` を持つ `NamedTuple`。`workspace` は現在の作業ディレクトリ。
"""
default_runtime_probe() = (
    platform = default_platform_probe(),
    kernel_release = default_kernel_probe(),
    workspace = pwd(),
)

"""
    workspace_check(runtime, workspace)

実行環境と作業場所を組み合わせて、教材の作業場所を確認する。

# 引数

- `runtime`: `classify_runtime` の分類結果。
- `workspace`: 観測した作業ディレクトリ。WSL2では `/home/` 配下を求める。

# 返り値

作業場所の判定・観測値・必要な対応を持つ `ObservedCheck`。ディレクトリの移動や作成はしない。
"""
function workspace_check(runtime, workspace)
    path = normpath(String(workspace))
    linux_home_workspace = startswith(path, "/home/") && length(path) > length("/home/")
    outside_linux_home = runtime.kind == :wsl2 && !linux_home_workspace
    passed = runtime.passed && !outside_linux_home
    observed = "pwd: $path"

    if outside_linux_home
        if path == "/mnt/c" || startswith(path, "/mnt/c/")
            observed *= " (Windows側/mnt/c)"
        else
            observed *= " (WSL2のLinux filesystem外)"
        end
        action = "WSL2では学生リポジトリを /home/<user>/... にcloneし、Linux filesystem内でF00を実行してください。"
    elseif !runtime.passed
        action = runtime.action
    else
        action = ""
    end
    ObservedCheck(:workspace, passed, observed, action)
end

"""
    default_command_probe(program, arguments)

PATH上のコマンドを実行し、利用可能かを観測する。

# 引数

- `program`: 実行するコマンド名。
- `arguments`: コマンドへ渡す引数の列。

# 返り値

`available` と `detail` を持つ `NamedTuple`。コマンドがない場合や終了コードが0以外の場合は `available=false`。標準出力・標準エラーを表示用の文字列にまとめる。
"""
function default_command_probe(program, arguments)
    executable = Sys.which(program)
    isnothing(executable) &&
        return (available = false, detail = "$(program)がPATH上に見つかりません")

    stdout = IOBuffer()
    stderr = IOBuffer()
    command = Cmd([executable, arguments...])
    process = run(pipeline(ignorestatus(command), stdout = stdout, stderr = stderr))
    output = strip(join(filter(!isempty, [String(take!(stdout)), String(take!(stderr))])))
    detail = isempty(output) ? "$(program)の終了コード: $(process.exitcode)" : output
    (available = process.exitcode == 0, detail = detail)
end

"""
    collect_preflight(;
        version_probe = () -> VERSION,
        command_probe = default_command_probe,
        runtime_probe = default_runtime_probe,
    )

実行環境・作業場所・Julia・Gitの観測結果をまとめる。

# 引数

- `version_probe`: Juliaの版を返す関数。既定は実行中の `VERSION`。
- `command_probe`: コマンドの `available` と `detail` を返す検査関数。
- `runtime_probe`: OS・カーネル名・作業場所の観測値を返す関数。

Juliaは要求する1.13系の版以上か検査する。検査関数で生じた例外は呼出し元へ伝わる。

# 返り値

4つの `ObservedCheck` を持つ `PreflightReport`。進捗は更新しない。
"""
function collect_preflight(;
    version_probe = () -> VERSION,
    command_probe = default_command_probe,
    runtime_probe = default_runtime_probe,
)
    runtime_input = runtime_probe()
    runtime_result = classify_runtime(runtime_input.platform, runtime_input.kernel_release)
    runtime_check = ObservedCheck(
        :runtime,
        runtime_result.passed,
        runtime_result.observed,
        runtime_result.action,
    )
    workspace = workspace_check(runtime_result, runtime_input.workspace)

    julia_version = version_probe()
    julia_check = ObservedCheck(
        :julia,
        julia_version >= REQUIRED_JULIA_VERSION &&
            (julia_version.major, julia_version.minor) ==
            (REQUIRED_JULIA_VERSION.major, REQUIRED_JULIA_VERSION.minor),
        string(julia_version),
        "JuliaupでJulia 1.13系をインストールして選択し、この確認を再実行してください。",
    )

    git_probe = command_probe("git", ["--version"])
    git_check = ObservedCheck(
        :git,
        git_probe.available,
        String(git_probe.detail),
        "Gitをインストールし、GitコマンドをPATHから実行できることを確認してください。",
    )

    PreflightReport(runtime_check, workspace, julia_check, git_check)
end

"""
    parse_preflight_arguments(arguments)

F00の手動確認オプションを解釈する。

# 引数

- `arguments`: `--confirm-vscode`、`--confirm-github`、`--confirm-agent <製品名>` の引数列。

重複・不明な引数・agentの値不足・未対応製品名は `ArgumentError` で拒否する。

# 返り値

`vscode_confirmed`、`github_confirmed`、`agent` を持つ `NamedTuple`。指定していない確認は `false`、agentは `nothing`。
"""
function parse_preflight_arguments(arguments)
    vscode_confirmed = false
    github_confirmed = false
    agent = nothing
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--confirm-vscode"
            vscode_confirmed &&
                throw(ArgumentError("--confirm-vscodeは1回だけ指定できます"))
            vscode_confirmed = true
            index += 1
        elseif argument == "--confirm-github"
            github_confirmed &&
                throw(ArgumentError("--confirm-githubは1回だけ指定できます"))
            github_confirmed = true
            index += 1
        elseif argument == "--confirm-agent"
            isnothing(agent) || throw(ArgumentError("--confirm-agentは1回だけ指定できます"))
            index == length(arguments) &&
                throw(ArgumentError("--confirm-agentには製品名が必要です"))
            candidate = arguments[index + 1]
            candidate in SUPPORTED_AGENTS || throw(
                ArgumentError(
                    "未対応のAIエージェント「$candidate」です。copilot、codex、amazon-qから選んでください",
                ),
            )
            agent = candidate
            index += 2
        else
            throw(ArgumentError("「preflight」の不明な引数です: $argument"))
        end
    end
    (; vscode_confirmed, github_confirmed, agent)
end

"""
    print_observed_check(io, label, check)

1つの観測結果と、必要な対応を表示する。

# 引数

- `io`: 表示先のIO。
- `label`: 項目名。
- `check`: 表示する `ObservedCheck`。

# 返り値

合格項目では `true`、要対応項目では `nothing`。`io` に判定と観測値を出力し、要対応なら対応方法も出力する。
"""
function print_observed_check(io, label, check)
    status = check.passed ? "PASS" : "NEEDS SETUP"
    println(io, "  [$status] $label: $(check.observed)")
    check.passed || println(io, "    対応: $(check.action)")
end

"""
    print_preflight(
        io,
        report;
        vscode_confirmed = false,
        github_confirmed = false,
        agent = nothing,
    )

機械観測と受講生の手動確認を分けて表示する。

# 引数

- `io`: 表示先のIO。
- `report`: `collect_preflight` で得た観測結果。
- `vscode_confirmed`: VS Codeの手動確認を済ませたか。
- `github_confirmed`: GitHubの手動確認を済ませたか。
- `agent`: 確認した対応AIエージェント名。未確認なら `nothing`。

# 返り値

`nothing`。観測項目と手動確認項目を `io` へ表示する。進捗の保存はしない。
"""
function print_preflight(
    io,
    report;
    vscode_confirmed = false,
    github_confirmed = false,
    agent = nothing,
)
    println(io, "端末で確認した項目")
    print_observed_check(io, "Runtime", report.runtime)
    print_observed_check(io, "作業ディレクトリ (pwd)", report.workspace)
    print_observed_check(io, "Julia", report.julia)
    print_observed_check(io, "Git", report.git)
    println(io)
    println(io, "手動確認")
    println(
        io,
        "  [$(vscode_confirmed ? "CONFIRMED" : "NOT CONFIRMED")] VS Code・Julia拡張機能 (WindowsはRemote - WSL)",
    )
    println(
        io,
        "  [$(github_confirmed ? "CONFIRMED" : "NOT CONFIRMED")] GitHubへのサインインとリポジトリへのアクセス",
    )
    agent_label = isnothing(agent) ? "未指定" : agent
    println(
        io,
        "  [$(isnothing(agent) ? "NOT CONFIRMED" : "CONFIRMED")] 正式対応AIエージェント: $agent_label",
    )
    nothing
end

"""
    run_f00_preflight(
        root;
        report = collect_preflight(),
        vscode_confirmed = false,
        github_confirmed = false,
        agent = nothing,
        persist_progress = save_progress,
        io = stdout,
    )

すべての観測と手動確認が完了した場合に、F00の進捗をF01へ進める。

# 引数

- `root`: 学生リポジトリのルート。
- `report`: 4項目の観測結果。既定では実際の環境を観測する。
- `vscode_confirmed`: VS Codeの手動確認結果。
- `github_confirmed`: GitHubの手動確認結果。
- `agent`: 確認した対応AIエージェント名。
- `persist_progress`: 進捗の保存関数。
- `io`: 診断と進捗判定の表示先。

未対応agentや、更新を許可しない進捗状態は `ArgumentError` で拒否する。未確認項目がある場合は進捗を変更しない。保存関数の失敗は呼出し元へ伝わる。

# 返り値

完了時は `true`、未確認項目が残る場合は `false`。初期F00の完了時だけ進捗ファイルを更新し、すでにF01へ進んだ状態は保持する。
"""
function run_f00_preflight(
    root;
    report = collect_preflight(),
    vscode_confirmed = false,
    github_confirmed = false,
    agent = nothing,
    persist_progress = save_progress,
    io = stdout,
)
    !isnothing(agent) &&
        !(agent in SUPPORTED_AGENTS) &&
        throw(ArgumentError("未対応のAIエージェント「$agent」です"))
    print_preflight(io, report; vscode_confirmed, github_confirmed, agent)

    observed_pass =
        report.runtime.passed &&
        report.workspace.passed &&
        report.julia.passed &&
        report.git.passed
    manual_pass = vscode_confirmed && github_confirmed && !isnothing(agent)
    if !(observed_pass && manual_pass)
        println(io)
        println(
            io,
            "F00の進捗は更新されませんでした。すべての機械観測と3項目の手動確認を完了してください。",
        )
        return false
    end

    progress_path = joinpath(root, "course_progress.toml")
    state = load_progress(progress_path)
    if state.current == "F01" && state.completed == ["F00"]
        println(io)
        println(io, "F00は完了済みです。現在の課題はF01のままです。")
        return true
    end
    state.current == "F00" && isempty(state.completed) ||
        throw(ArgumentError("F00の事前診断は初期F00進捗だけを更新できます"))

    advanced = ProgressState(state.schema_version, state.ordered, ["F00"], "F01")
    persist_progress(progress_path, advanced)
    println(io)
    println(
        io,
        "F00が完了しました。現在の課題はF01です。F00用のbranchやPRは作成しないでください。",
    )
    true
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    F00Environment.print_preflight(stdout, F00Environment.collect_preflight())
    println()
    println(
        "このスクリプトは診断専用です。https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/F00.html に従い、scripts/course.jlからF00を完了してください。",
    )
end
