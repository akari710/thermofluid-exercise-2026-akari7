# N08/N09: Neumann発展の提供済みテスト。境界・適合条件・平均0を独立に確認します。

using Test
include("neumann.jl");
include("extension_support.jl")
const N = NeumannExtension
@testset "混合境界から純Neumann" begin
    errors = Float64[]
    for (nx, ny) in OFFICIAL_GRIDS, id in NEUMANN_CASES
        p = neumann_problem(id, nx, ny)
        result = N.solve_neumann(p.u0, p.f, p.dx, p.dy, p.bc)
        @test result.converged &&
              extension_residual(result.u, p.f, p.dx, p.dy, p.bc) <=
              result.threshold + 3e-11
        exact = copy(p.exact)
        if is_pure(p.bc)
            exact.-=weighted_mean(exact, p.dx, p.dy)
            @test abs(weighted_mean(result.u, p.dx, p.dy)) < 1e-13
            pin = result.u .- result.u[1, 1]
            @test extension_residual(pin, p.f, p.dx, p.dy, p.bc) <= result.threshold + 3e-11
            @test all(
                maximum(abs, a - b) < 2e-13 for (a, b) in zip(
                    values(normal_derivatives(pin, p.dx, p.dy)),
                    values(normal_derivatives(result.u, p.dx, p.dy)),
                )
            )
        end
        if id == "mixed_poisson_zero_flux"
            @test result.u[div(nx, 2), 1] > 0.5 &&
                  result.u[1, 1] == 0 &&
                  result.u[end, end] == 4
            @test maximum(abs, result.u - exact) <= 2result.threshold + 1e-10
        elseif id == "pure_poisson_zero_flux"
            push!(errors, sqrt(sum(abs2, result.u - exact) / (nx * ny)))
        else
            @test maximum(abs, result.u - exact) < 0.025
        end
    end
    @test all(v->1.8 <= v <= 2.2, log2.(errors[1:2] ./ errors[2:3]))
    p = neumann_problem("pure_poisson_flux", 17, 13)
    f = p.f .+ 1
    copyf = copy(f)
    @test_throws ArgumentError N.solve_neumann(p.u0, f, p.dx, p.dy, p.bc)
    @test f == copyf
end
