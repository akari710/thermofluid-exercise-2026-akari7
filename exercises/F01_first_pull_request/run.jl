# F01の挨拶プログラム。入力検証 → 学生が実装する文字列生成 → 表示の順に読みます。
# 実装対象は student_greeting のTODOです。main と直接実行の入口は提供済みです。

module F01FirstPullRequest

export main, student_greeting

"""
    student_greeting(name::AbstractString)::String

名前の前後の空白を除き、挨拶文字列を作る。

# 引数

- `name`: 挨拶に使う名前。空白を除いた結果が空なら `ArgumentError`。

# 返り値

実装後は `"Hello, <name>!"` の形式の `String`。

# 受講生のToDo

空白を除いた名前で挨拶文字列を返す。配布時は `未実装 F01: student_greeting` で停止する。
"""
function student_greeting(name::AbstractString)::String
    normalized_name = strip(name)
    isempty(normalized_name) && throw(ArgumentError("名前を空にはできません"))

    # TODO(F01): `Hello, <normalized_name>!`を返す処理を実装する。
    error("未実装 F01: student_greeting")
end

"""
    main(name::AbstractString = "student"; io = stdout)

挨拶を生成し、指定した表示先へ出力する。

# 引数

- `name`: 挨拶する名前。既定は `"student"`。
- `io`: 表示先のIO。既定は標準出力。

実行には `student_greeting` のTODO実装が必要。名前の検証や未実装エラーは呼出し元へ伝わる。

# 返り値

生成した挨拶の `String`。同じ文字列を `io` へ改行付きで表示する。
"""
function main(name::AbstractString = "student"; io = stdout)
    message = student_greeting(name)
    println(io, message)
    message
end

end


if abspath(PROGRAM_FILE) == @__FILE__
    name = isempty(ARGS) ? "student" : only(ARGS)
    F01FirstPullRequest.main(name)
end
