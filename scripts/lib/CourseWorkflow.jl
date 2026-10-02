module CourseWorkflow

using TOML

export ORDERED_UNITS, TASK_IDS_BY_UNIT, ProgressState, load_progress, save_progress,
       tests_to_run, validate_transition, UNIT_DIRECTORIES, unit_directory, units_to_test, require_unit_assets

const ORDERED_UNITS = [
    "F00", "F01", "F02", "F03-F04",
    "N01", "N02", "N03", "N04", "N05-N06", "N07", "N08-N09",
]
const TASK_IDS_BY_UNIT = Dict(unit => [unit] for unit in ORDERED_UNITS)
TASK_IDS_BY_UNIT["F03-F04"] = ["F03", "F04"]
TASK_IDS_BY_UNIT["N05-N06"] = ["N05", "N06"]
TASK_IDS_BY_UNIT["N08-N09"] = ["N08", "N09"]

const PROGRESS_KEYS = Set(["schema_version", "ordered", "completed", "current"])

struct ProgressState
    schema_version::Int
    ordered::Vector{String}
    completed::Vector{String}
    current::String
end

function _string_vector(value, field)
    value isa AbstractVector || throw(ArgumentError("`$field`は提出単位の配列である必要があります"))
    all(item -> item isa AbstractString, value) ||
        throw(ArgumentError("`$field`には文字列の提出単位だけを指定してください"))
    String[String(item) for item in value]
end

function _validate(state::ProgressState)
    state.schema_version == 2 ||
        throw(ArgumentError("未対応の`schema_version`です: $(state.schema_version)"))
    state.ordered == ORDERED_UNITS ||
        throw(ArgumentError("`ordered`には全提出単位を授業順で1回ずつ指定してください"))
    length(unique(state.completed)) == length(state.completed) ||
        throw(ArgumentError("`completed`に重複した提出単位があります"))
    state.current in state.ordered ||
        throw(ArgumentError("`current`に不明な提出単位があります: $(state.current)"))
    state.current in state.completed &&
        throw(ArgumentError("`current`の提出単位を`completed`にも含めることはできません"))

    current_index = findfirst(==(state.current), state.ordered)
    expected_completed = state.ordered[1:(current_index - 1)]
    state.completed == expected_completed ||
        throw(ArgumentError("`completed`には`current`より前の提出単位を順番どおり指定してください"))
    nothing
end

function load_progress(path)
    values = TOML.parsefile(path)
    Set(keys(values)) == PROGRESS_KEYS ||
        throw(ArgumentError("進捗ファイルには`schema_version`、`ordered`、`completed`、`current`だけを含めてください"))

    schema_version = values["schema_version"]
    schema_version isa Integer && !(schema_version isa Bool) ||
        throw(ArgumentError("`schema_version`は整数である必要があります"))
    current = values["current"]
    current isa AbstractString || throw(ArgumentError("`current`は文字列の提出単位である必要があります"))

    state = ProgressState(
        Int(schema_version),
        _string_vector(values["ordered"], "ordered"),
        _string_vector(values["completed"], "completed"),
        String(current),
    )
    _validate(state)
    state
end

function save_progress(path, state::ProgressState)
    _validate(state)
    directory = dirname(abspath(path))
    temporary_path, output = mktemp(directory; cleanup=false)
    try
        TOML.print(output, Dict(
            "schema_version" => state.schema_version,
            "ordered" => state.ordered,
            "completed" => state.completed,
            "current" => state.current,
        ); sorted=true)
        close(output)
        mv(temporary_path, path; force=true)
    finally
        isopen(output) && close(output)
        rm(temporary_path; force=true)
    end
    nothing
end

const UNIT_DIRECTORIES = Dict(
    "F00" => "F00_environment",
    "F01" => "F01_first_pull_request",
    "F02" => "F02_julia_arrays_and_tests",
    "F03-F04" => "F03-F04_vector_calculus",
    "N01" => "N01_linear_advection",
    "N02" => "N02_nonlinear_advection",
    "N03" => "N03_diffusion",
    "N04" => "N04_advection_diffusion",
    "N05-N06" => "N05-N06_common_package_2d_advection",
    "N07" => "N07_2d_advection_diffusion",
    "N08-N09" => "N08-N09_laplace_poisson",
)

unit_directory(id) = joinpath("exercises", UNIT_DIRECTORIES[id])

function require_unit_assets(root, id)
    missing = String[]
    directory = unit_directory(id)
    required = ["run.jl", "tests.jl", "learning_log.md"]
    id == "F03-F04" && push!(required, "F03.jl")
    id in ("N01", "N02", "N03", "N04") && push!(required, "provided_support.jl")
    if id == "N05-N06"
        append!(required,["N05.jl","simulate.jl","analyze.jl","plot.jl","provided_support.jl"])
        isfile(joinpath(root,"src","N06Advection.jl")) || push!(missing,"src/N06Advection.jl")
    end
    if id == "N07"
        append!(required,["simulate.jl","analyze.jl","plot.jl","provided_support.jl"])
        isfile(joinpath(root,"src","N07Transport.jl")) || push!(missing,"src/N07Transport.jl")
    end
    if id == "N08-N09"
        append!(required,["simulate.jl","analyze.jl","plot.jl","provided_support.jl"])
        isfile(joinpath(root,"src","N08N09Elliptic.jl")) || push!(missing,"src/N08N09Elliptic.jl")
    end
    for name in required
        relative = joinpath(directory, name)
        isfile(joinpath(root, relative)) || push!(missing, relative)
    end
    isempty(missing) || throw(ArgumentError(
        "$id の教材が揃っていません: $(join(missing, ", "))。現在の課題を続け、教員に配布状況を確認してください",
    ))
    nothing
end

function units_to_test(state::ProgressState)
    _validate(state)
    filter(!=("F00"), vcat(state.completed, [state.current]))
end

function tests_to_run(state::ProgressState)
    _validate(state)
    units = vcat(state.completed, [state.current])
    reduce(vcat, (TASK_IDS_BY_UNIT[unit] for unit in units); init=String[])
end

function validate_transition(state::ProgressState, id)
    _validate(state)
    id in state.ordered ||
        throw(ArgumentError("不明な提出単位です: `$id`"))
    id != "F00" ||
        throw(ArgumentError("F00の環境確認にはpreflightを使用してください"))
    nothing
end

end
