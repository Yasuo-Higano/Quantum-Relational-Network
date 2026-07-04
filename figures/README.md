# 論文図 (v9.3)

`make_figures.py` が results/ の機械可読 JSON (一次ソース) から 5 点の SVG を生成する。
**図中の全数値は一次ソースと照合され、1 つでも不一致なら [FAIL] で exit 1 する** —
この照合が v9.2 の「世代ラベルの綱渡り」を発見した装置でもある (docs/uft-v9.2.md §1)。

## 再現手順

```bash
python3 -m venv figenv && figenv/bin/pip install matplotlib numpy
figenv/bin/python figures/make_figures.py            # SVG を figures/ に生成
figenv/bin/python figures/make_figures.py --png DIR  # 検収用 PNG も DIR に出力
```

(CLAUDE.md 追加規則「人間向け資料は python 可」による。シミュレーション本体は Rust。)

## 図の一覧

| ファイル | 論文 | 内容 | 一次ソース |
|---|---|---|---|
| `v62_landscape.svg` | anomaly 図 1 | 理論空間スペクトル (v6.2 で手書き SVG) | results/v62_atlas.txt |
| `fig_controls_map.svg` | anomaly 図 2 | 陰性対照の地図: 最小性は {カイラル性, 全因子帯電, Witten, SU(3)³} が、一意性は U(1)³ が担う | results/v62_atlas.json |
| `fig_u1_staircase.svg` | anomaly 図 3 | U(1) の階段: Y → B−L → 3 本目なし (対照 355) | results/v62_atlas.json, v71_twou1.json, v82_threeu1.json |
| `fig_zeromode_wilson.svg` | yukawa 図 1 | ゼロモード局在と Wilson 線の剛体平行移動 | 再計算 (results/v72_geomfn.txt と照合) |
| `fig_attainable_ratios.svg` | yukawa 図 2 | 到達可能な質量比集合: 単一 T² の床 3×10⁻³ と T²×T² の 2 乗抑制 | 再計算 (v72/v92 と照合) |
| `fig_geometry_lnz.svg` | yukawa 図 3 | 幾何選択の証拠 (質量のみ/全 9 量 × 両ラベル規約) + 緊張の解消 | results/v92_labelstab.json, v81/v91/v72 JSON |

図ラベルは投稿用に英語。世代ラベルは v9.2 の安定規約。
