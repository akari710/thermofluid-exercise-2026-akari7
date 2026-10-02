include(joinpath(@__DIR__, "lib", "CourseWorkflow.jl"))
include(joinpath(@__DIR__, "lib", "ResultLimits.jl"))
include(joinpath(@__DIR__, "lib", "CourseTests.jl"))
include(joinpath(@__DIR__, "..", "exercises", "F00_environment", "run.jl"))

using .CourseWorkflow
using .F00Environment
using .ResultLimits
using .CourseTests

const SLUGS = Dict(
    "F01" => "first-pull-request",
    "F02" => "julia-arrays-and-tests",
    "F03-F04" => "vector-calculus-numerical-differentiation",
    "N01" => "linear-advection",
    "N02" => "nonlinear-advection",
    "N03" => "diffusion",
    "N04" => "advection-diffusion",
    "N05-N06" => "common-package-2d-advection",
    "N07" => "2d-advection-diffusion",
    "N08-N09" => "laplace-poisson",
)

const USAGE = """
使い方:
  julia --project=. $(joinpath("scripts", "course.jl")) preflight [--confirm-vscode --confirm-github --confirm-agent <copilot|codex|amazon-q>]
  julia --project=. $(joinpath("scripts", "course.jl")) start <ID>
  julia --project=. $(joinpath("scripts", "course.jl")) status
  julia --project=. $(joinpath("scripts", "course.jl")) check-results

例:
  julia --project=. $(joinpath("scripts", "course.jl")) start F02
  julia --project=. $(joinpath("scripts", "course.jl")) start F03-F04

必要な教材が配布済みの任意の提出単位を開始できます。F00の環境確認にはpreflightを使用してください。
"""

git_output(root, arguments...) = readchomp(Cmd(`git $(arguments)`; dir=root))

function require_preflight(root)
    branch = git_output(root, "branch", "--show-current")
    branch == "main" ||
        throw(ArgumentError("courseコマンドはmainから開始してください。現在のbranch: $branch"))
    isempty(git_output(root, "status", "--porcelain")) ||
        throw(ArgumentError("作業ツリーに未commitの変更があります。git status --shortで確認してください"))
    nothing
end

function require_result_limits(root)
    violations = check_result_limits(root)
    isempty(violations) || throw(ArgumentError(join(violations, '\n')))
    nothing
end

function require_local_tests(root)
    run_course_tests(root; policy=:advance) || error("共通検査が失敗しました。進捗は更新しません")
    nothing
end

function show_status(root)
    state = load_progress(joinpath(root, "course_progress.toml"))
    completed = isempty(state.completed) ? "なし" : join(state.completed, ", ")
    println("現在の提出単位: $(state.current)")
    println("通過済み: $completed")
    if state.current == "F00"
        println("環境診断: julia --project=. scripts/course.jl preflight")
        return
    end
    for id in TASK_IDS_BY_UNIT[state.current]
        println("課題ページ: https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/$id.html")
    end
    directory = unit_directory(state.current)
    println("開くフォルダ: $directory/")
    if state.current == "N05-N06"
        println("編集するファイル: src/、analyze.jl、plot.jl、tests.jl、learning_log.md")
    elseif state.current == "N07"
        println("編集するファイル: src/N07Transport.jl、analyze.jl、tests.jl、learning_log.md（表示はplot.jl）")
    elseif state.current == "N08-N09"
        println("編集するファイル: src/N08N09Elliptic.jl、tests.jl、learning_log.md（計算simulate.jl、解析analyze.jl、表示plot.jl）")
        println("授業内: julia --project=. $(joinpath(directory,"tests.jl")) N08-check（N08の4関数のみ。全体完了ではありません）")
        println("N08単独: julia --project=. $(joinpath(directory,"run.jl")) N08")
        println("全体: julia --project=. $(joinpath(directory,"run.jl")) all")
        println("発展はextensions/の独立コマンドです。必修の完了には要求しません")
        println("src編集後: julia --project=. exercises/N05-N06_common_package_2d_advection/N05.jl verify → Pkg.test()")
    else
        println("編集するファイル: run.jl、tests.jl、learning_log.md")
    end
    state.current in ("N05-N06", "N07", "N08-N09") && println("共通コード: src/")
    println("実行: julia --project=. $(joinpath(directory, "run.jl"))")
    println("テスト: julia --project=. -e 'using Pkg; Pkg.test()'")
end

function start_exercise(root, id; persist_progress=save_progress)
    progress_path = joinpath(root, "course_progress.toml")
    state = load_progress(progress_path)
    require_preflight(root)
    validate_transition(state, id)
    require_unit_assets(root, id)
    require_local_tests(root)

    slug = get(SLUGS, id, nothing)
    isnothing(slug) &&
        throw(ArgumentError("課題ID $id のbranch slugが設定されていません"))
    branch = "exercise/$id-$slug"
    run(Cmd(`git switch -c $branch`; dir=root))
    advanced = ProgressState(
        state.schema_version,
        state.ordered,
        state.ordered[1:(findfirst(==(id), state.ordered) - 1)],
        id,
    )
    try
        persist_progress(progress_path, advanced)
    catch persistence_error
        try
            run(Cmd(`git switch main`; dir=root))
            run(Cmd(`git branch -D $branch`; dir=root))
        catch rollback_error
            throw(CompositeException([persistence_error, rollback_error]))
        end
        rethrow()
    end

    println("$(id)をbranch $(branch)で開始しました。")
    show_status(root)
    println("作業をcommitしたら、次のコマンドでpushしてください:")
    println("  git push -u origin $branch")
    println("その後、ホスティングサービスで$(branch)のpull requestを作成してください。")
end

function main(
    arguments=ARGS;
    root=pwd(),
    preflight_report=nothing,
    preflight_collector=collect_preflight,
    persist_progress=save_progress,
    io=stdout,
)
    if arguments == ["--help"] || arguments == ["-h"] || isempty(arguments)
        print(USAGE)
        return 0
    end

    command = first(arguments)
    if command == "preflight"
        confirmations = parse_preflight_arguments(arguments[2:end])
        report = isnothing(preflight_report) ? preflight_collector() : preflight_report
        completed = run_f00_preflight(
            root;
            report,
            persist_progress,
            io,
            confirmations...,
        )
        attempted_completion = length(arguments) > 1
        return attempted_completion && !completed ? 1 : 0
    elseif command == "status" && length(arguments) == 1
        show_status(root)
    elseif command == "check-results" && length(arguments) == 1
        require_result_limits(root)
        println("結果ファイルのサイズ上限を確認しました。")
    elseif command == "start" && length(arguments) == 2
        start_exercise(root, arguments[2])
    else
        throw(ArgumentError("無効なコマンドです\n$USAGE"))
    end
    0
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(main())
    catch exception
        println(stderr, "エラー: ", sprint(showerror, exception))
        exit(1)
    end
end
