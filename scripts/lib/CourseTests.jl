module CourseTests

using Test
using TOML
using ..CourseWorkflow
using ..ResultLimits

export run_course_tests

# The child writes a result only after reaching and finishing the exercise boundary.
# Startup, package loading, interruption and abnormal exit must not become warnings.
function exercise_worker(path, report)
    kind = "assertion"
    detail = ""
    status = "passed"
    try
        @testset "$(basename(dirname(path))) / tests.jl" begin
            try
                Base.include(Main, path)
            catch exception
                cause = exception
                while cause isa LoadError
                    cause = cause.error
                end
                cause isa InterruptException && throw(cause)
                kind = cause isa Base.Meta.ParseError ? "syntax" : "exception"
                detail = sprint(showerror, exception, catch_backtrace())
                rethrow()
            end
        end
    catch exception
        exception isa Test.TestSetException || rethrow()
        status = "failed"
        isempty(detail) && (detail = sprint(showerror, exception))
    end
    open(report, "w") do io
        TOML.print(io, Dict("status" => status, "kind" => kind, "detail" => detail))
    end
    nothing
end

function julia_process(root, code, arguments...)
    Cmd(Cmd([Base.julia_cmd().exec..., "--startup-file=no", "--project=$root",
             "-e", code, arguments...]); dir=root)
end

function require_success(command, label)
    process = run(ignorestatus(command))
    success(process) || error("$label が失敗しました（実行基盤または共通検査）。進捗は更新しません")
end

function collect_result(root, worker, path, label)
    mktempdir() do temporary
        report = joinpath(temporary, "result.toml")
        require_success(julia_process(root, worker, root, path, report), "$label runner")
        isfile(report) || error("$label runnerの結果を回収できません。進捗は更新しません")
        result = TOML.parsefile(report)
        Set(keys(result)) == Set(["status", "kind", "detail"]) &&
            result["status"] in ("passed", "failed") &&
            result["kind"] in ("assertion", "syntax", "exception") &&
            result["detail"] isa String || error("$label runnerの結果が不正です")
        result
    end
end

function run_course_tests(root; policy=:current, io=stderr)
    policy in (:current, :advance) || throw(ArgumentError("不明なテスト方針: $policy"))
    root = abspath(root)
    state = load_progress(joinpath(root, "course_progress.toml"))
    units = units_to_test(state)
    required = ["Project.toml", "Manifest.toml", "src/ThermofluidExercise.jl", "src/N06Advection.jl",
                "src/N07Transport.jl", "scripts/lib/CourseWorkflow.jl",
                "scripts/lib/ResultLimits.jl", "scripts/lib/CourseTests.jl",
                "exercises/F00_environment/run.jl", "test/f00_preflight_test.jl",
                "test/course_transition_test.jl"]
    append!(required, [joinpath(unit_directory(unit), "tests.jl") for unit in units])
    for relative in required
        isfile(joinpath(root, relative)) || error("必須ファイルがありません: $relative")
    end
    for unit in units
        require_unit_assets(root, unit)
    end
    require_limits() = begin
        violations = check_result_limits(root)
        isempty(violations) || error(join(violations, '\n'))
    end
    require_limits()
    worker = """
        using ThermofluidExercise
        include(joinpath(ARGS[1], "scripts", "lib", "CourseWorkflow.jl"))
        include(joinpath(ARGS[1], "scripts", "lib", "ResultLimits.jl"))
        include(joinpath(ARGS[1], "scripts", "lib", "CourseTests.jl"))
        using .CourseWorkflow
        CourseTests.exercise_worker(ARGS[2], ARGS[3])
        """
    common = collect_result(root, "using HDF5, Plots\n" * worker,
                            joinpath(root, "test", "f00_preflight_test.jl"), "F00共通検査")
    all_passed = common["status"] == "passed"
    if !all_passed
        println(io, "エラー: F00共通検査が失敗しました。残りの検査を続行します。")
        println(io, common["detail"])
    end
    for unit in units
        path = joinpath(root, unit_directory(unit), "tests.jl")
        result = collect_result(root, worker, path, unit)
        if result["status"] == "failed"
            allowed = policy == :advance || unit != state.current
            level = allowed ? "警告" : "エラー"
            println(io, "$level: $unit / $(relpath(path, root)) ($(result["kind"])) が失敗しました。残りの検査を続行します。")
            println(io, result["detail"])
            allowed || (all_passed = false)
        end
    end
    require_limits()
    all_passed
end

end
