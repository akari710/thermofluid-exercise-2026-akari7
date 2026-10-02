using Test

@testset "任意の課題への直接開始" begin
    for current_index in eachindex(ORDERED_UNITS), target in ORDERED_UNITS[2:end]
        state = ProgressState(2, ORDERED_UNITS, ORDERED_UNITS[1:current_index-1], ORDERED_UNITS[current_index])
        @test validate_transition(state, target) === nothing
    end
    initial = ProgressState(2, ORDERED_UNITS, String[], "F00")
    @test_throws ArgumentError validate_transition(initial, "unknown")
    @test_throws ArgumentError validate_transition(initial, "F00")
    mktempdir() do root
        for (index, target) in enumerate(ORDERED_UNITS[2:end])
            path = joinpath(root, "course_progress.toml")
            save_progress(path, ProgressState(2, ORDERED_UNITS, ORDERED_UNITS[1:index], target))
            restored = load_progress(path)
            @test restored.current == target
            @test restored.completed == ORDERED_UNITS[1:index]
        end
    end
end
