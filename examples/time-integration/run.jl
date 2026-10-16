module TimeIntegrationExample
using LinearAlgebra, TOML
include("provided_methods.jl")

"""
    rk2_step(f, y, t, dt)

f(y,t)を用いる陽的中点法の1ステップ。
状態は有限実数または非空のVector{Float64}、時刻は有限、dtは有限正。
配列入力を保持し、独立した配列を返す。
配布済みrhs(f,y,t)は右辺の形と有限性を確認する。
"""
function rk2_step(f,y,t,dt)
    validate_step(y,t,dt)
    # TODO(N03): 陽的中点法の二段更新を実装する（今年度は任意比較）。
    error("未実装 N03: rk2_step (任意比較)")
end

include("provided_support.jl")
end

if abspath(PROGRAM_FILE) == @__FILE__
    arguments = filter(!=("--cn"),ARGS)
    length(arguments)<=1 || throw(ArgumentError("使い方: run.jl [出力先] [--cn]"))
    directory = isempty(arguments) ? TimeIntegrationExample.default_output_dir() : only(arguments)
    TimeIntegrationExample.main(;output_dir=directory,include_cn="--cn" in ARGS)
end
