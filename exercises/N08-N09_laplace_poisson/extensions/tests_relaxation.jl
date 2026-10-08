# N08/N09: 緩和法発展の提供済みテスト。逐次更新と境界保持を確認します。

using Test, ThermofluidExercise
include("relaxation.jl");
include("extension_support.jl")
const R = RelaxationExtension
@testset "GSとSOR" begin
    p = problem("N09", 17, 13)
    a = copy(p.u0)
    b = copy(a)
    R.gauss_seidel_step!(a, p.f, p.dx, p.dy)
    R.sor_step!(b, p.f, p.dx, p.dy; omega = 1.)
    @test a == b
    @test_throws ArgumentError R.sor_step!(b, p.f, p.dx, p.dy; omega = 2.)
    counts = Int[]
    for method in (:gauss_seidel, :sor)
        result = R.solve_relaxation(p.u0, p.g, p.f, p.dx, p.dy; method)
        @test result.converged &&
              independent_residual(result.u, p.f, p.dx, p.dy) <= result.threshold + 2e-11
        push!(counts, result.iterations)
    end
    @test counts[2] < counts[1]
end
