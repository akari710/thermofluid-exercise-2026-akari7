# N07: 提供済みの独立した期待値による必須確認 → 受講生が編集する自作テスト。

module N07Tests
using Test, ThermofluidExercise
const N = ThermofluidExercise.N07
include("provided_support.jl")
include("analyze.jl")
# 配布済みの必須テスト。自作欄と区別し、入力・期待値を保持する。
@testset "N07 independent faces, half-cell walls and heat" begin
    a = Float64[1 2 4 3; 5 7 6 8; 9 10 12 11]
    saved = copy(a)
    b = similar(a)
    bc = N.channel_boundaries(:dirichlet)
    f = N.flux_buffers(a)
    @test N.thermal_stable_timestep(1., 0., 0.05, 0.5, 0.25, bc)≈0.8 / 5.0
    @test N.thermal_fluxes!(f, a, 0.5, 0.25; cx=1., cy=0., kappa=0.05, bc)===f
    @test f.adv_x == [1 1 1 1; 1 2 4 3; 5 7 6 8; 9 10 12 11]
    @test f.diff_x≈[0 -0.2 -0.6 -0.4; -0.4 -0.5 -0.2 -0.5; -0.4 -0.3 -0.6 -0.3; 0 0 0 0]
    @test f.adv_y == zeros(3, 5)
    @test f.diff_y≈[-0.4 -0.2 -0.4 0.2 1.2; -2 -0.4 0.2 -0.4 3.2; -3.6 -0.2 -0.4 0.2 4.4]
    r = N.thermal_step!(b, a, 0.01, 0.5, 0.25; cx=1., cy=0., kappa=0.05, bc)
    @test r.temperature===b
    @test collect(r.boundary_rates.advective) == [1., -10.5, 0., 0.]
    @test collect(r.boundary_rates.diffusive)≈[-0.3, 0., -3., -4.4]
    @test b≈[1.0 1.994 3.908 2.922; 4.856 6.872 5.992 7.752; 8.776 9.942 11.844 10.766] atol =
        1e-14
    @test 0.5 * 0.25 * sum(b - a)≈0.01 * (-17.2) atol = 1e-14
    @test a == saved
    for wall in (:dirichlet, :insulated)
        bc = N.channel_boundaries(wall)
        q0 = 0.125sum(a)
        total = 0.
        u = copy(a)
        v = similar(u)
        for _ in 1:20
            r = N.thermal_step!(v, u, 0.01, 0.5, 0.25; cx=1., cy=0., kappa=0.05, bc)
            rates = r.boundary_rates
            total+=0.01 * (sum(rates.advective) + sum(rates.diffusive))
            u, v = v, u
        end
        @test 0.125sum(u) - q0≈total atol = 1e-12
        wall == :insulated &&
            @test r.boundary_rates.diffusive.south == r.boundary_rates.diffusive.north == 0.
    end
    for kind in (:periodic, :insulated, :dirichlet)
        bc = N.boundaries(kind; value=kind == :dirichlet ? 2. : 0.)
        N.thermal_step!(
            b,
            fill(2., 3, 4),
            0.01,
            0.5,
            0.25;
            cx=0.,
            cy=0.,
            kappa=0.05,
            bc,
        )
        @test b≈fill(2., 3, 4)
    end
    bc = N.boundaries(:dirichlet)
    bound = 1 / (0.05 * (3 / 0.5 ^ 2 + 3 / 0.25 ^ 2))
    impulse = zeros(3, 4)
    impulse[1, 1] = 1.
    N.thermal_step!(b, impulse, bound, 0.5, 0.25; cx=0., cy=0., kappa=0.05, bc)
    @test minimum(b) >= -1e-15
    fill!(b, -99.)
    @test_throws ArgumentError N.thermal_step!(
        b,
        a,
        bound * (1 + 1e-8),
        0.5,
        0.25;
        cx=0.,
        cy=0.,
        kappa=0.05,
        bc,
    )
    @test all(==(-99.), b)
    @test_throws ArgumentError N.thermal_step!(
        b,
        a,
        1 / (0.1 * (1 / 0.5 ^ 2 + 1 / 0.25 ^ 2)),
        0.5,
        0.25;
        cx=0.,
        cy=0.,
        kappa=0.05,
        bc,
    )
    # Signed periodic flow at both ends, old-state use and heat cancellation.
    bc = N.boundaries(:periodic)
    for (cx, cy, kappa) in ((-1., 0.5, 0.), (0., 0., 0.05), (0., 0., 0.))
        r = N.thermal_step!(b, a, 0.01, 0.5, 0.25; cx, cy, kappa, bc)
        @test sum(b)≈sum(a) atol = 1e-12
        @test sum(r.boundary_rates.advective) + sum(r.boundary_rates.diffusive)≈0. atol =
            1e-13
        if kappa == 0 && cx != 0
            @test b[3, 1]≈8.88
        elseif cx == cy == kappa == 0
            @test b == a
        end
    end
end
@testset "N07 two-component signed advective Burgers" begin
    u = [1. -2. 0. 0.5; -0.5 2. -1. 0.; 3. -1. 0.2 -2.]
    v = [-0.5 1. 2. -1.; 2. -1. 0. 0.5; 0. 0.5 -2. 1.]
    originals = (copy(u), copy(v))
    un = similar(u)
    vn = similar(v)
    dt = 0.001
    dx = 0.5
    dy = 0.25
    nu = 0.05
    r = N.burgers_step!(un, vn, u, v, dt, dx, dy; nu)
    @test r.u===un && r.v===vn
    @test un≈[
        0.9953 -1.967 -0.01736 0.4989;
        -0.4891 1.9662 -0.99396 -0.0035;
        2.9707 -0.98944 0.17892 -1.97454
    ] atol = 1e-14
    @test vn≈[
        -0.4946 0.9851 1.9876 -0.9925;
        1.9815 -0.9841 -0.0044 0.4996;
        0.0135 0.4974 -1.9696 0.9763
    ] atol = 1e-14
    @test u == originals[1] && v == originals[2]
    @test N.burgers_stable_timestep(u, v, nu, dx, dy)≈0.8 / (9.0 + 2.0)
    @test N.burgers_stable_timestep(2u, 2v, nu, dx, dy) <
        N.burgers_stable_timestep(u, v, nu, dx, dy)
    a = fill(2., 3, 4)
    b = fill(-1., 3, 4)
    N.burgers_step!(un, vn, a, b, 0.001, dx, dy; nu=0.)
    @test un == a && vn == b
    N.burgers_step!(un, vn, u, zeros(3, 4), dt, dx, dy; nu=0.)
    @test all(iszero, vn)
    arrays = [copy(u) for _ in 1:4]
    for (i, j) in ((1, 3), (1, 2))
        args = copy(arrays)
        args[j] = args[i]
        @test_throws ArgumentError N.burgers_step!(args..., dt, dx, dy; nu)
    end
    fill!(un, -99.)
    fill!(vn, -88.)
    @test_throws ArgumentError N.burgers_step!(un, vn, u, v, 1., dx, dy; nu)
    @test all(==(-99.), un) && all(==(-88.), vn)
    @test_throws ArgumentError N.burgers_stable_timestep(u, v, -1., dx, dy)
end

@testset "N07 配布済み熱収支解析" begin
    b = N07Analysis.heat_budget(
        [0., 0.25, 1.],
        [2., 4., 5.],
        [3. 4.; -1. -2.; 0. 0.; 0. 0.],
        [0.5 0.2; 0. 0.; -0.2 -0.4; -0.3 -0.8],
    )
    @test b.net_input≈[2., 1.]
    @test b.cumulative_input≈[0., 2., 3.]
    @test b.residual≈zeros(3) atol = 1e-14
    @test_throws ArgumentError N07Analysis.heat_budget(
        [0., 0.],
        [1., 1.],
        zeros(4, 1),
        zeros(4, 1),
    )
end
@testset "N07 自作: 開境界2条件の辺別熱収支" begin
    # TODO(自作): 複数ステップで辺別・移流/拡散別の積算、熱量変化、断熱壁を検証する。
    # 入力と独立期待値、許容値の根拠を記す。公式出力は上書きしない。
    initial = [2.0 4.0 1.0;
        3.0 6.0 5.0;
        7.0 2.0 8.0;
        4.0 9.0 3.0]

    nx, ny = size(initial)
    dx, dy = 0.4, 0.2
    cx, cy = 0.7, 0.0
    kappa = 0.04
    dt = 0.01
    nsteps = 12
    inlet_temperature = 1.0
    wall_temperature = 0.0

    initial_heat = dx * dy * sum(initial)

    # 12ステップの加算と，熱量の差を取る際の丸め誤差を許容する．
    tolerance = 1e-12 * max(1.0, abs(initial_heat))

    for wall in (:dirichlet, :insulated)
        @testset "$wall" begin
            bc = N.channel_boundaries(wall)
            old = copy(initial)
            new = similar(old)

            expected_adv = zeros(4)
            expected_diff = zeros(4)
            actual_adv = zeros(4)
            actual_diff = zeros(4)

            for _ in 1:nsteps
                # 更新前の温度から，西・東・南・北の流入レートを求める．
                west_adv = cx * dy * ny * inlet_temperature
                east_adv = -cx * dy * sum(old[end, :])
                west_diff = kappa * dy / (dx / 2) *
                    (ny * inlet_temperature - sum(old[1, :]))

                south_diff = 0.0
                north_diff = 0.0
                if wall == :dirichlet
                    wall_factor = kappa * dx / (dy / 2)
                    south_diff = wall_factor *
                        (nx * wall_temperature - sum(old[:, 1]))
                    north_diff = wall_factor *
                        (nx * wall_temperature - sum(old[:, end]))
                end

                # レートに dt を掛けて，輸送された熱量を積算する．
                expected_adv .+= dt .* [west_adv, east_adv, 0.0, 0.0]
                expected_diff .+= dt .* [west_diff, 0.0, south_diff, north_diff]

                result = N.thermal_step!(
                    new, old, dt, dx, dy; cx, cy, kappa, bc,
                )
                rates = result.boundary_rates
                actual_adv .+= dt .* collect(rates.advective)
                actual_diff .+= dt .* collect(rates.diffusive)

                if wall == :insulated
                    @test iszero(rates.diffusive.south)
                    @test iszero(rates.diffusive.north)
                end

                old, new = new, old
            end

            for edge in 1:4
                @test actual_adv[edge] ≈ expected_adv[edge] atol=tolerance rtol=0
                @test actual_diff[edge] ≈ expected_diff[edge] atol=tolerance rtol=0
            end

            heat_change = dx * dy * sum(old) - initial_heat
            expected_input = sum(expected_adv) + sum(expected_diff)
            @test heat_change ≈ expected_input atol=tolerance rtol=0
        end
    end
end
@testset "N07 自作: 温度とBurgersの解析解収束" begin
    # TODO(自作): 周期温度の拡散・合成と二成分Burgersを3格子で計算し、誤差と両区間の次数を検証する。
    grids = (24, 48, 96)
    final_time = 0.2
    kappa = 0.005
    nu = 0.005
    bc = N.boundaries(:periodic)

    rms_error(numerical, exact) =
        sqrt(sum(abs2, numerical .- exact) / length(exact))

    function temperature_error(n, cx, cy)
        h = 1.0 / n
        centers = ((1:n) .- 0.5) .* h

        function exact_temperature(t)
            amplitude = 0.3 * exp(-8 * pi^2 * kappa * t)
            return [
                1.0 + amplitude * sinpi(2 * (x + y - (cx + cy) * t))
                for x in centers, y in centers
            ]
        end

        old = exact_temperature(0.0)
        new = similar(old)

        dt_max = N.thermal_stable_timestep(cx, cy, kappa, h, h, bc)
        nsteps = ceil(Int, final_time / dt_max)
        dt = final_time / nsteps

        for _ in 1:nsteps
            N.thermal_step!(new, old, dt, h, h; cx, cy, kappa, bc)
            old, new = new, old
        end

        return rms_error(old, exact_temperature(final_time))
    end

    function burgers_error(n)
        h = 1.0 / n
        centers = ((1:n) .- 0.5) .* h
        u_mean, v_mean = 0.65, -0.35

        function exact_wave(t)
            amplitude = 0.15 * exp(-8 * pi^2 * nu * t)
            return [
                amplitude * sinpi(2 * (x + y - (u_mean + v_mean) * t))
                for x in centers, y in centers
            ]
        end

        initial_wave = exact_wave(0.0)
        u_old = u_mean .+ initial_wave
        v_old = v_mean .- initial_wave
        u_new = similar(u_old)
        v_new = similar(v_old)

        t = 0.0
        while t < final_time
            dt_max = N.burgers_stable_timestep(u_old, v_old, nu, h, h)
            dt = min(dt_max, final_time - t)

            N.burgers_step!(u_new, v_new, u_old, v_old, dt, h, h; nu)
            u_old, u_new = u_new, u_old
            v_old, v_new = v_new, v_old
            t += dt
        end

        final_wave = exact_wave(final_time)
        u_error = rms_error(u_old, u_mean .+ final_wave)
        v_error = rms_error(v_old, v_mean .- final_wave)
        return u_error, v_error
    end

    function check_convergence(errors, p_min, p_max)
        @test all(isfinite, errors)
        @test errors[1] > errors[2] > errors[3] > 0.0

        # 格子幅は各区間で半分になるので，分母は log(2)．
        for i in 1:2
            order = log2(errors[i] / errors[i+1])
            @test p_min <= order <= p_max
        end
    end

    @testset "周期温度: 拡散" begin
        errors = [temperature_error(n, 0.0, 0.0) for n in grids]
        check_convergence(errors, 1.8, 2.2)
    end

    @testset "周期温度: 移流拡散" begin
        errors = [temperature_error(n, 0.65, -0.35) for n in grids]
        check_convergence(errors, 0.8, 1.2)
    end

    @testset "周期Burgers: u成分とv成分" begin
        errors = [burgers_error(n) for n in grids]
        check_convergence(first.(errors), 0.8, 1.2)
        check_convergence(last.(errors), 0.8, 1.2)
    end
end
@testset "N07 保存結果と来歴" begin
    @test check_complete()
end
end
