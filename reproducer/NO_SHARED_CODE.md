# クリーンルーム規約 — 何を「独立外部再現」と数えるか (v27.4 凍結)

本リポジトリの最大の不足は内部計算量ではなく、**外部に固定された再現単位と
独立の検証**である (PROMPT/9 の結論)。その不足を正確に埋めるため、
「独立外部再現」の定義を機械可読に固定する ([replications.yml](../replications.yml))。

## 数えるための 6 条件 (全て必須)

```yaml
different_author: true            # 本リポジトリの作者と別人 (別 AI を含む — 下記)
independent_repository: true      # 別リポジトリで公開されている
no_shared_numerical_kernel: true  # 数値カーネル (固有値分解・積分器・列挙器) を共有しない
protocol_frozen_before_run: true  # 走行前に SPEC/TOLERANCES のコミットを固定した
commit_hash_recorded: true        # 再現実装のコミットハッシュが報告に記録されている
result_including_failures_public: true  # 失敗・不一致を含む全結果が公開されている
```

## 数えないもの (algorithmic diversity — 価値はあるが独立外部再現ではない)

- **同一作者**が Python・Julia・Rust などで再実装したもの (本リポジトリ内の
  numpy/LAPACK 照合・Prolog 独立推論・Lean 形式化はこの類 —
  claims.yml の independence: algorithmically_diverse / same_author_clean_room)。
- **同一 AI (または同系列のモデル)** が仕様と実装の両方を書いたもの。本リポジトリ
  の実装は AI 支援で書かれているため、同じ AI による「独立実装」は独立性の
  実質を欠く。
- 本リポジトリの数値カーネルを流用・移植したもの (バグも一緒に移植される)。
- 走行後に許容誤差・プロトコルを調整したもの (覗き見)。

## 独立性が壊れやすい点 (再現者への注意)

- 単位 A の正準化 (並べ替え・共役・符号・gcd) は**自分で設計**すること。本
  リポジトリの正準形 SHA-256 と合わなくても、多重項集合が一致すれば再現成立。
- 単位 B は閉形式の評価だけでなく、有限鎖の直接対角化 (別経路) を推奨する。
- 疑問点は issue で質問してよい — ただし回答は仕様の明確化に限り、
  実装コードは提供しない。

## 現在の計数

**独立外部再現 = 0** (この数字は v271_core_audit [S4/S5] と v274_reproducer が
機械監査しており、成立時には監査の期待値ごと更新される)。
