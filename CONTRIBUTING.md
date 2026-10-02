# コントリビューション

## 対応環境

Windowsで課題を実行する場合は，WSL2 Ubuntu 24.04 のLinux側でJulia，Git，SSH，VS Code Remote - WSL，Julia拡張機能，エージェントを使います．
リポジトリはLinux側の `/home/<user>/...` に置きます．
macOSとnative Linuxは各OS側の環境を使います．

## 履修者の課題branchとPR

`main`から課題branchを作り，変更・テスト・学習ログをそろえてから`main`へのpull request（PR）を作成します．
PRではテンプレートに沿って，diffとGitHub Actionsを確認してください．

教材更新は`upstream`から更新用branchへ通常のmergeで取り込み，更新PRを**Create a merge commit**で統合します．
現在課題のテスト失敗は修正し，過去課題の警告も課題IDと原因を確認してください．
作業環境を再準備した場合は，F00確認後に進捗をcommit・ローカルmergeしてから，`start TASK_ID`で再開する課題を直接指定します．
詳細は[課題ワークフロー](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/workflow.html#previous-exercises)と[配布更新](https://t2lab-it.github.io/thermofluid-exercise-2026/guides/commands.html#material-updates)を参照してください．

## 外部からの提案

外部の提案は，このリポジトリを`fork`してbranchを作り，forkからPRを作成してください．
提案には，次を記載します．

- 再利用した素材の出典URL
- 採用・改変した内容と変更点
- 実行したテスト，整形確認，または数値的な検証の根拠

## 再利用とライセンス

コード，設定，テスト，スクリプトへの寄与はMIT License，READMEの文章，公開学習ログ，学生作成図への寄与はCC BY 4.0で提供することに同意したうえでPRを送ってください．
第三者素材を含める場合は，その利用条件，出典，必要な表示を確認してください．

## 公開してはいけない内容

秘密情報，個人情報，成績情報，LMSの提出内容，非公開URL，アクセス情報，または他者を傷つける内容をPRやIssueに書かないでください．
理解度チェックのAI対話全文は，UTF-8のテキストファイルとしてLETUSへ提出しますが，リポジトリ，commit，PR，学習ログには含めません．
合理的配慮に関する情報，その他のAI利用の全文ログ，生の会話記録は，リポジトリに含めません．
学習ログにはAI利用の要約，または「利用なし」のみを記録してください．
個人情報や有害な内容に関する相談は，LMSの個別連絡または授業で案内された教員連絡先から非公開の連絡をしてください．
