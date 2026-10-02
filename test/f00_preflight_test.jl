using Test

const F00_REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(F00_REPO_ROOT, "exercises", "F00_environment", "run.jl"))
using .F00Environment
using .CourseWorkflow

include(joinpath(@__DIR__, "course_transition_test.jl"))

const WSL2_KERNEL = "5.15.167.4-microsoft-standard-WSL2"

function f00_report(;
    platform=:linux,
    kernel_release="6.8.0-generic",
    workspace="/home/student/thermofluid-exercise",
    julia_version=v"1.13.0",
    git_ok=true,
)
    command_probe = function (program, arguments)
        (program, arguments) == ("git", ["--version"]) ||
            error("unexpected command probe: $program $(join(arguments, ' '))")
        (available=git_ok, detail=git_ok ? "git version test" : "git not found")
    end
    collect_preflight(;
        version_probe=() -> julia_version,
        command_probe,
        runtime_probe=() -> (; platform, kernel_release, workspace),
    )
end

function initial_f00_repo()
    root = mktempdir()
    save_progress(
        joinpath(root, "course_progress.toml"),
        ProgressState(2, ORDERED_UNITS, String[], "F00"),
    )
    root
end

@testset "F00 environment preflight" begin
    @testset "supported runtimes and WSL2 workspace" begin
        for (platform, kernel_release, kind, passed) in (
            (:linux, WSL2_KERNEL, :wsl2, true),
            (:linux, "6.8.0-generic", :native_linux, true),
            (:macos, "Darwin", :macos, true),
            (:windows, "", :native_windows, false),
            (:linux, "4.4.0-Microsoft", :wsl1, false),
            (:linux, "", :unknown_linux, false),
        )
            report = f00_report(; platform, kernel_release)
            @test (classify_runtime(platform, kernel_release).kind,
                   report.runtime.passed, report.workspace.passed) == (kind, passed, passed)
        end
        for workspace in ("/mnt/c/course", "/tmp/course")
            report = f00_report(; kernel_release=WSL2_KERNEL, workspace)
            @test !report.workspace.passed
            @test occursin(workspace, report.workspace.observed)
            @test occursin("/home/<user>", report.workspace.action)
        end
    end

    @testset "Julia version and Git diagnostics" begin
        report = f00_report()
        @test report.julia.passed && report.julia.observed == "1.13.0"
        @test report.git.passed

        for julia_version in (v"1.13.0", v"1.13.1", v"1.13.99")
            compatible = f00_report(; julia_version)
            @test compatible.julia.passed
            @test compatible.julia.observed == string(julia_version)
        end
        for julia_version in (v"1.12.9", v"1.13.0-rc1", v"1.14.0-rc1", v"1.14.0", v"2.0.0")
            wrong_version = f00_report(; julia_version)
            @test !wrong_version.julia.passed
            @test occursin("1.13", wrong_version.julia.action)
        end
        missing_git = f00_report(git_ok=false)
        @test !missing_git.git.passed
        @test occursin("Git", missing_git.git.action)
    end

    @testset "manual confirmation arguments" begin
        @test parse_preflight_arguments(String[]) ==
            (vscode_confirmed=false, github_confirmed=false, agent=nothing)
        for agent in ("copilot", "codex", "amazon-q")
            @test parse_preflight_arguments([
                "--confirm-agent", agent, "--confirm-github", "--confirm-vscode",
            ]) == (vscode_confirmed=true, github_confirmed=true, agent=agent)
        end
        for arguments in (
            ["--confirm-agent", "claude"],
            ["--confirm-agent"],
            ["--confirm-agent", "copilot", "--confirm-agent", "codex"],
            ["--unknown"],
        )
            @test_throws ArgumentError parse_preflight_arguments(arguments)
        end
    end

    @testset "each gate is required before changing progress" begin
        root = initial_f00_repo()
        progress_path = joinpath(root, "course_progress.toml")
        before = read(progress_path, String)
        report = f00_report()

        for (vscode_confirmed, github_confirmed, agent) in (
            (false, true, "codex"),
            (true, false, "codex"),
            (true, true, nothing),
        )
            @test !run_f00_preflight(
                root; report, vscode_confirmed, github_confirmed, agent, io=IOBuffer(),
            )
            @test read(progress_path, String) == before
        end
        for failed_report in (
            f00_report(julia_version=v"1.14.0"),
            f00_report(git_ok=false),
            f00_report(platform=:windows),
            f00_report(kernel_release=WSL2_KERNEL, workspace="/mnt/c/course"),
        )
            @test !run_f00_preflight(
                root; report=failed_report, vscode_confirmed=true,
                github_confirmed=true, agent="codex", io=IOBuffer(),
            )
            @test read(progress_path, String) == before
        end

        output = IOBuffer()
        @test run_f00_preflight(
            root; report, vscode_confirmed=true, github_confirmed=true,
            agent="codex", io=output,
        )
        state = load_progress(progress_path)
        @test (state.completed, state.current) == (["F00"], "F01")
        @test occursin("F00が完了しました。現在の課題はF01です。", String(take!(output)))
    end

    @testset "persistence failure and repeated completion preserve progress" begin
        root = initial_f00_repo()
        progress_path = joinpath(root, "course_progress.toml")
        before = read(progress_path, String)
        report = f00_report()
        @test_throws ErrorException run_f00_preflight(
            root; report, vscode_confirmed=true, github_confirmed=true, agent="copilot",
            persist_progress=(path, state) -> error("injected persistence failure"),
            io=IOBuffer(),
        )
        @test read(progress_path, String) == before
        @test run_f00_preflight(
            root; report, vscode_confirmed=true, github_confirmed=true,
            agent="copilot", io=IOBuffer(),
        )
        after = read(progress_path, String)
        @test run_f00_preflight(
            root; report, vscode_confirmed=true, github_confirmed=true, agent="amazon-q",
            persist_progress=(path, state) -> error("completed progress must not be rewritten"),
            io=IOBuffer(),
        )
        @test read(progress_path, String) == after
        @test !isdir(joinpath(root, ".git"))
    end
end

@testset "F00 course CLI rejects missing Git without changing progress" begin
    root = initial_f00_repo()
    progress_path = joinpath(root, "course_progress.toml")
    before = read(progress_path, String)
    # This process checks CLI wiring and exit status; it does not need numerical JIT optimization.
    command = Cmd(Cmd([
        Base.julia_cmd().exec...,
        "--startup-file=no", "--compile=min", "-O0",
        "--project=$F00_REPO_ROOT",
        joinpath(F00_REPO_ROOT, "scripts", "course.jl"),
        "preflight", "--confirm-vscode", "--confirm-github", "--confirm-agent", "codex",
    ]); dir=root)
    output = IOBuffer()
    process = run(pipeline(addenv(ignorestatus(command), "PATH" => mktempdir()), stdout=output))
    text = String(take!(output))
    @test process.exitcode == 1
    @test occursin("[NEEDS SETUP] Git", text)
    @test occursin("F00の進捗は更新されませんでした", text)
    @test read(progress_path, String) == before
end
