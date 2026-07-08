//! v11.4 動的テンソル分解の第一級化 — 「どの分解が自然か」は状態が決める
//!
//! 残高 7 (core 移行) の最終項目。v1.3/v2.0 の追補「テンソル分解は任意ではない」を、
//! core の読み出し `readout_natural_partition` (サイト間相互情報量の貪欲最大マッチング)
//! として第一級化する。同じ読み出しが、状態ごとに**異なる**自然な分解を発見すること —
//! つまり分解が力学的・状態依存であること — を 3 つの core 模型で検証する:
//!
//!   [1] RingChain 基底状態      → 空間的に隣接する対 (局所性の再発見; v6.4 と整合)
//!   [2] TfdPair (熱場二重)      → L–R の鏡像対 (空間ではなく「橋」が自然な分解)
//!   [3] GrowingChain 到着直後   → 到着ボンド対そのもの (成長の単位が分解の単位)
//!       同・時間発展後          → 到着対の分解は溶けて組み変わる (分解は動的)
//!
//! [1][2][3] は厳密な構造の検査 (全対一致)。[3'] は分解が発展で変わることの実演。
//! これで core 移行地図 (RingChain / TfdPair / GrowingChain / PacketRing + 分解読み出し)
//! は残高なく閉じる。

use uft_sim::*;

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

fn main() {
    self_test();
    println!("=== v11.4 動的テンソル分解の第一級化: 分解は状態が決める ===\n");

    // ---- [1] RingChain: 自然な分解 = 空間隣接対 ----
    println!("[1] RingChain (n=42) 基底状態 — 自然な分解は空間的な隣接対か");
    let ring = RingChain { n: 42 };
    let s_ring = ring.init();
    let (pairs_r, frac_r) = s_ring.readout_natural_partition();
    let n = 42usize;
    let adjacent = pairs_r
        .iter()
        .all(|&(i, j)| (i + 1) % n == j || (j + 1) % n == i);
    println!(
        "    対 21 組: 全て隣接 (|i−j| = 1 mod n)  {}  (捕捉 MI 割合 {:.3})",
        pass(adjacent),
        frac_r
    );

    // ---- [2] TfdPair: 自然な分解 = L–R 鏡像対 ----
    println!("\n[2] TfdPair (各鎖 n=21, β=1) — 自然な分解は L–R の鏡像対か");
    let tfd = TfdPair { n: 21, beta: 1.0 };
    let s_tfd = tfd.init();
    let (pairs_t, frac_t) = s_tfd.readout_natural_partition();
    let mirror = pairs_t
        .iter()
        .all(|&(i, j)| (i < 21 && j == i + 21) || (j < 21 && i == j + 21));
    println!(
        "    対 21 組: 全て鏡像 (x_L ↔ x_R)  {}  (捕捉 MI 割合 {:.3})",
        pass(mirror),
        frac_t
    );
    println!("    — 空間の隣接ではなく「橋」が分解の単位: ER=EPR の分解版");

    // ---- [3] GrowingChain: 到着直後は到着対、発展後は組み変わる ----
    println!("\n[3] GrowingChain (n0=6 → 12; 到着 3 回) — 分解は動的か");
    let gc = GrowingChain { n_max: 12 };
    let mut s = gc.init(6);
    for m in 0..3usize {
        s = gc.arrive_pair_vacuum(&s, 6 + 2 * m);
    }
    let (pairs_g, frac_g) = s.readout_natural_partition();
    // 到着した 3 対 (6,7), (8,9), (10,11) がそのまま対として発見されるか
    let arrived_ok = [(6usize, 7usize), (8, 9), (10, 11)].iter().all(|&(a, b)| {
        pairs_g
            .iter()
            .any(|&(i, j)| (i, j) == (a, b) || (j, i) == (a, b))
    });
    println!(
        "    到着直後: 到着ボンド対 (6,7)(8,9)(10,11) が分解に現れる  {}  (捕捉 {:.3})",
        pass(arrived_ok),
        frac_g
    );
    // 発展させると到着対の純度が下がり、分解が組み変わる
    let s2 = gc.evolve_active(&s, 12, 30.0);
    let (pairs_g2, frac_g2) = s2.readout_natural_partition();
    let arrived_after = [(6usize, 7usize), (8, 9), (10, 11)]
        .iter()
        .filter(|&&(a, b)| {
            pairs_g2
                .iter()
                .any(|&(i, j)| (i, j) == (a, b) || (j, i) == (a, b))
        })
        .count();
    println!(
        "    発展後 (t=30): 到着対のうち分解に残るのは {}/3 — 対が溶けて組み変わった (捕捉 {:.3})",
        arrived_after, frac_g2
    );
    let dynamic = arrived_after < 3;
    println!("    分解は動的 (発展で変わる)  {}", pass(dynamic));

    // ---- [4] 健全性: 3 状態とも純粋なガウス状態のまま ----
    println!("\n[4] 健全性 (core の共通検査)");
    let pd = s_ring
        .purity_defect()
        .max(s_tfd.purity_defect())
        .max(s.purity_defect());
    println!(
        "    purity_defect 最大 {:.2e} (< 1e-9)  {}",
        pd,
        pass(pd < 1e-9)
    );

    let all_ok = adjacent && mirror && arrived_ok && dynamic && pd < 1e-9;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-CORE-004".into())),
        ("ring_adjacent".into(), Json::Bool(adjacent)),
        ("ring_captured".into(), Json::Num(frac_r)),
        ("tfd_mirror".into(), Json::Bool(mirror)),
        ("tfd_captured".into(), Json::Num(frac_t)),
        ("growing_arrival_pairs".into(), Json::Bool(arrived_ok)),
        (
            "growing_pairs_surviving_evolution".into(),
            Json::Int(arrived_after as i64),
        ),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v114_tensorcore.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {}", pass(all_ok));
    println!("\n結論: 同じ読み出しが、円環では空間の隣接を、熱場二重では橋 (鏡像対) を、");
    println!("      成長宇宙では到着の単位を「自然な分解」として発見し、発展はそれを");
    println!("      組み替える。テンソル分解は入力ではなく読み出しである — 残高 7 は閉じた。");
    if !all_ok {
        std::process::exit(1);
    }
}
