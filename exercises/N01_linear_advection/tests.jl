# N01のテスト枠。必須テストと自作テストの両方を受講生が実装します。
# 各TODOで入力と独立した期待値を選び、配列全体・境界・安定性を確認します。

using Test

if !isdefined(Main, :N01LinearAdvection)
    include(joinpath(@__DIR__, "run.jl"))
end

using .N01LinearAdvection

@testset "N01 必須テスト" begin
    # TODO(必須): 座標、base、plateau、区間端を自分で選び、rectangular_initial_conditionが返す配列全体を手計算した期待値と比較する。
    # 1. 矩形状初期分布
    x = [0.0, 0.5, 1.0, 1.5, 2.0]

    u = rectangular_initial_condition(
        x;
        base=1.0,
        plateau=3.0,
        plateau_start=0.5,
        plateau_end=1.5
    )

    @test u == [1.0, 3.0, 3.0, 3.0, 1.0]

    # TODO(必須): 小さなu_old、移流速度、時間刻み、格子幅を選び、upwind_step!後のu_new配列全体を手計算した期待値と比較する。
    # 2. 風上差分法
    u_old = [1.0, 2.0, 4.0, 3.0, 1.0]
    u_new = similar(u_old)

    upwind_step!(u_new, u_old, 1.0, 0.5, 1.0)

    @test u_new == [1.0, 1.5, 3.0, 3.5, 1.0]

    # TODO(必須): 同じ入力でcentered_step!後のu_new配列全体を手計算した期待値と比較する。風上差分との違いが現れる入力を選ぶ。
    # 3. 中心差分法
    u_new = similar(u_old)

    centered_step!(u_new, u_old, 1.0, 0.5, 1.0)

    @test u_new == [1.0, 1.25, 3.75, 3.75, 1.0]

    # TODO(必須): apply_boundary!が左端と右端を課題の境界条件どおりに書き換えることを、配列全体の期待値で確かめる。
    # 4. 境界条件
    u = [10.0, 2.0, 3.0, 4.0, 20.0]

    apply_boundary!(u; left_value=1.0)

    @test u == [1.0, 2.0, 3.0, 4.0, 4.0]

    # TODO(必須): simulateで安定な風上差分と不安定な中心差分を同じ条件で計算し、最終値の範囲または振幅から両者の違いを確かめる。条件と判定基準は自分で書く。
    # 5. 安定性の比較
    upwind = simulate(
        scheme=:upwind,
        nx=81,
        c=1.0,
        cfl=0.5,
        t_final=0.5
    )

    centered = simulate(
        scheme=:centered,
        nx=81,
        c=1.0,
        cfl=0.5,
        t_final=0.5
    )

    @test minimum(upwind.u) >= 1.0 - 1e-10 &&
        maximum(upwind.u) <= 2.0 + 1e-10 &&
        (minimum(centered.u) < 1.0 - 1e-6 ||
            maximum(centered.u) > 2.0 + 1e-6)

end

@testset "N01 自作テスト" begin
    # TODO(自作): 必須とは異なる初期条件、CFL、格子数、またはAPIの性質を一つ選び、期待する結果を自分で書く。PNG、TOML、ファイルサイズは対象にしない。
    # 矩形区間内に格子点がない場合
    x = [0.0, 1.0, 2.0, 3.0]

    u = rectangular_initial_condition(
        x;
        base=1.0,
        plateau=3.0,
        plateau_start=1.2,
        plateau_end=1.8
    )

    @test u == [1.0, 1.0, 1.0, 1.0]
end
