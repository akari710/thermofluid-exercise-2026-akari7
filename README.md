# 熱流体力学演習（2026）学生用リポジトリ

数式・実装・テストを対応させ，数値結果と限界を説明するための個人課題リポジトリです．

<!-- contract-section: student_created_repository -->
## 初回の準備

- Julia 1.13系
- Git
- VS Code
- GitHub Copilot，OpenAI Codex，Amazon Q Developerのいずれか一つ

## 対応環境

Windows では，演習に必要な Julia，Git，SSH，VS Code，Julia 拡張機能，コーディングエージェントを WSL2 の Ubuntu 24.04 側で実行します．
WSL2 の学生リポジトリは `/home/<user>/...` に SSH で複製し，Windows の VS Code は Remote - WSL でそのLinux側フォルダーを開いてください．
macOS と native Linux は，それぞれのOS側の Julia，Git，SSH，VS Code，エージェントを使います．

GitHubで自分をOwnerとする空の公開リポジトリ `thermofluid-exercise-2026-<自分のusername>` を作成します．
README・`.gitignore`・ライセンスの初期生成のチェックを外して作成します．
[学生用配布リポジトリ](https://github.com/t2lab-it/thermofluid-exercise-student-2026)をSSHで通常cloneして履歴を保持し，配布元を`upstream`，個人の提出先を`origin`にします．
以下の`OWNER/REPOSITORY`と`YOUR_COURSE_REPOSITORY`を自分の値に置き換えます．

```fish
git clone git@github.com:t2lab-it/thermofluid-exercise-student-2026.git YOUR_COURSE_REPOSITORY
cd YOUR_COURSE_REPOSITORY
git remote rename origin upstream
git remote add origin git@github.com:OWNER/REPOSITORY.git
git push -u origin main
```

別PCで個人リポジトリをcloneしたら，`upstream`を再登録します．
[環境診断](https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/F00.html)と[Git・GitHubの準備](https://t2lab-it.github.io/thermofluid-exercise-2026/setup/git-github.html)に沿って進めてください．

```fish
julia --project=. scripts/course.jl preflight --confirm-vscode --confirm-github --confirm-agent codex
```

## 課題を開く

[公開課題ページ](https://t2lab-it.github.io/thermofluid-exercise-2026/)から，対応するローカルの課題フォルダを開きます．
通常編集するのは `run.jl`（実装），`tests.jl`（自分の確認），`learning_log.md`（記録）です．
実装箇所は`TODO(課題ID)`で探し，`error("未実装 課題ID: 関数識別子")`を自分の計算処理へ置き換えます．
課題ページに従って`src/`，`analyze.jl`，発展の`extensions/`も確認し，N06の分散の例外と選択・発展の範囲は[共通ガイド](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/workflow.html#implementation-todos)を参照してください．
課題ごとの `tests.jl` には，必須の数値・入出力テストと自分で設計するテストの TODO があります．

理解度チェックを行った授業では，AIとの対話全文を `understanding-check-<lesson-id>.txt` というUTF-8テキストに保存し，対応するLETUSの課題へ提出します．
課題内の `learning_log.md` には提出済みの確認と学んだ点を記録します．
公開する内容は[コントリビューション](CONTRIBUTING.md#公開してはいけない内容)を確認してください．

| 提出単位 | フォルダ | 内容 |
|---|---|---|
| F00 | `exercises/F00_environment/` | 環境診断 |
| F01 | `exercises/F01_first_pull_request/` | 最初のPR |
| F02 | `exercises/F02_julia_arrays_and_tests/` | 配列・関数・テスト |
| F03-F04 | `exercises/F03-F04_vector_calculus/` | ベクトル解析・数値微分 |
| N01 | `exercises/N01_linear_advection/` | 一次元線形移流 |
| N02 | `exercises/N02_nonlinear_advection/` | 一次元非線形移流 |
| N03 | `exercises/N03_diffusion/` | 一次元熱拡散 |
| N04 | `exercises/N04_advection_diffusion/` | 選択モデルの一次元移流拡散 |

```fish
julia --project=. scripts/course.jl status
```

F01の手動branch作成，F02以降の課題開始，提出順は[課題ワークフロー](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/workflow.html)を参照してください．

## 実行とテスト

リポジトリのルートで，課題ページに記載された `run.jl` を実行します．

```fish
julia --project=. -e 'using Pkg; Pkg.test()'
```

過去課題の失敗は警告として現在課題まで検査を続け，現在課題の失敗はエラーとして扱います．
詳しくは[テスト結果の読み方](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/testing.html#test-results)を参照してください．
`start TASK_ID`では課題の失敗を警告にして，指定した任意の提出単位を開始できます．
例として，`current = "F00"`から`julia --project=. scripts/course.jl start F02`でF02を直接開始できます．
F00の環境確認には`preflight`を使用します．
選択した課題より前を通過済みとして記録し，その課題以降は進捗記録から外します．
`completed`と通過済みは進捗の記録です．
課題の完了は，各課題ページの完了条件で確認してください．
結果がある課題では，実行時に課題内の `results/` が作られます．
図・数値を確認し，ログから参照してcommitします．

詳しい[コマンド一覧](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/commands.html)と[最終プロジェクトへの移行](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/final-project-handoff.html)も参照してください．

## 最終プロジェクト

通常課題の成果を保ったまま，学生自身が公開リポジトリを新規作成します．
AIと協働して，必要なコード，Julia環境，テスト，CI，再現手順を整えます．
実施形態（2人1組か1人1テーマか）は未定です．
リポジトリの作成単位と共同作業者の設定は，実施形態の決定後に案内します．
[最終プロジェクトの環境構築とコード移行](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/final-project-handoff.html)に従い，必要な課題コードとテストだけを移し，出典と参照元commitを記録してください．
