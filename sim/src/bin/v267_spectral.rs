//! v26.7 (II/II) v267_spectral — 有限サイズ spectral measure と massless pole 判定 (経路 B)
//!
//! 事前登録: paper/grav-vacuum-polarization-spec.md §7 (157ca53)。主成果物は broadened
//! plot ではなく**有限サイズ spectral measure** ρ_A^(N)(s) = Σ_n Z_n δ(s − s_n)
//! (s = ΔE², 粒子–正孔対の離散和)。pole 判定 (凍結): massless pole の主張は
//! s_n → 0 **かつ** Z_n → Z_* > 0 (residue の体積 scaling) を要する —
//! 「分母が零に近い」は証拠にならない。**自由場の登録予想 (v26.5): pole なし
//! (連続体のみ)** — 誘導重力の pole は 1/Π の構造から来るべきもの。
//!
//! 設定: q = 2π/16 ŷ を固定 (物理運動量固定の体積列 N ∈ {16, 32, 64}, j = N/16)、
//! m ∈ {0, 0.5}。チャネルは proofs/Projector.lean の ŷ 辞書:
//! D = (T_xx−T_zz)/√2 (spin-2 plus), S = (T_xx+T_zz)/√2 (P0s), L = T_yy (縦=純ゲージ)。
//! 対ごとの weight: Z_A = |⟨μ|T_A(q)|ν⟩|², s = (E_μ − E_ν)²。
//!
//! 検査 (凍結):
//!  [S0] 内部恒等: 対和 Σ 2Z_A/ΔE = ブロック Lehmann (v26.4–26.6 経路の χ_A) — 1e-12
//!       (対ごとの分解が χ 器械と同じ物理を刻んでいること)
//!  [S1] **f-sum rule (独立モーメント照合)**: Σ_n 2 ΔE_n Z_n^A =
//!       ⟨0|[T_A(q)†, [H, T_A(q)]]|0⟩ (二重交換子は 1 体演算子 — 相関行列の trace で
//!       独立に計算)。dense (N=8) と block (全 N) の両方で相対 1e-10
//!  [S2] 正値性: min Z_A ≥ 0 (構成的 — 実頂点の 2 乗) と s_min > 0
//!  [S3] threshold (m=0.5): s_min ≥ 4m² (厳密不等式 E(k) ≥ m の帰結 — 機械精度)
//!  [S4] massless edge (m=0): s_min(N)/q² の記録 — 連続体は光円錐 ω = q から始まる
//!       (branch 記録: 端が q² に近づくか)
//!  [S5] **pole 判定 (branch)**: W₁(N) = (最低クラスタの Σ Z)/V の体積 scaling。
//!       branch α: W₁ が N とともに減少 → 自由場に pole なし (登録予想の確認) /
//!       branch β: W₁ → Z_* > 0 → pole 様 (予想が外れる — それ自体を公表)。
//!       ゲートは「分解能」(順序が単調に確定すること) のみ — どちらの branch も判定 a
//!  [S6] 変異: 折返しスワップ落とし → S0/S1 が検出 (> 1e-4)
//!
//! 事前登録分岐: (a) S0–S3 PASS → spectral 器械認証、S5 の branch が主結果 /
//!   (b) S0/S1 FAIL → 対分解の誤り (dense が真) / (c) S5 分解能なし → 体積列の拡張。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

// ---------------- ブロック機構 (v26.6 の認証済み実装を写経) ----------------

fn block_h(n: usize, m: f64, cky: f64, ckz: f64) -> Vec<f64> {
    let dim = 4 * n;
    let mut h = vec![0.0f64; dim * dim];
    let id = |x: usize, c: usize| x + n * c;
    let add = |h: &mut Vec<f64>, a: usize, b: usize, t: f64| {
        h[b + a * dim] += t;
        h[a + b * dim] += t;
    };
    let ysgn = [1.0, -1.0, 1.0, -1.0];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        for c in 0..4 {
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            add(&mut h, id(x, c), id((x + 1) % n, c), 0.5 * tw);
            h[id(x, c) + id(x, c) * dim] += px * ysgn[c] * cky;
            if c == 0 || c == 2 {
                add(&mut h, id(x, c), id(x, c + 1), px * zsgn[c] * ckz);
            }
        }
        add(&mut h, id(x, 0), id(x, 3), px * m);
        add(&mut h, id(x, 1), id(x, 2), px * m);
    }
    h
}

/// q∥ŷ の頂点 (v26.4–26.6 と同一): which 1 = T_xx, 2 = T_yy, 3 = T_zz
fn vertex_qy(n: usize, ky: f64, ckz: f64, q: f64, sw: bool, which: usize) -> Vec<f64> {
    let dim = 4 * n;
    let id = |x: usize, c: usize| x + n * c;
    let tc = |c: usize| -> usize {
        if sw {
            [1usize, 0, 3, 2][c]
        } else {
            c
        }
    };
    let mut v = vec![0.0f64; dim * dim];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        for c in 0..4 {
            let kyc = ky + if c == 1 || c == 3 { PI } else { 0.0 };
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            if which == 1 {
                v[id((x + 1) % n, tc(c)) + id(x, c) * dim] += 0.5 * tw;
                v[id(x, tc(c)) + id((x + 1) % n, c) * dim] += 0.5 * tw;
            }
            if which == 2 {
                v[id(x, tc(c)) + id(x, c) * dim] += px * (kyc + q / 2.0).cos();
            }
            if which == 3 && (c == 0 || c == 2) {
                let coef = px * zsgn[c] * ckz;
                v[id(x, tc(c + 1)) + id(x, c) * dim] += coef;
                v[id(x, tc(c)) + id(x, c + 1) * dim] += coef;
            }
        }
    }
    v
}

/// チャネル頂点 D, S, L (実, B1→B2 の長方 dim×dim 格納 v[標的 + 元·dim])
fn channel_vertices(n: usize, ky: f64, ckz: f64, q: f64, sw: bool) -> [Vec<f64>; 3] {
    let ox = vertex_qy(n, ky, ckz, q, sw, 1);
    let oy = vertex_qy(n, ky, ckz, q, sw, 2);
    let oz = vertex_qy(n, ky, ckz, q, sw, 3);
    let dim = 4 * n;
    let r2 = (2.0f64).sqrt();
    let mut od = vec![0.0f64; dim * dim];
    let mut os = vec![0.0f64; dim * dim];
    for k in 0..dim * dim {
        od[k] = (ox[k] - oz[k]) / r2;
        os[k] = (ox[k] + oz[k]) / r2;
    }
    [od, os, oy]
}

/// 1 ブロック対の spectral 収集 + f-sum RHS。
/// 戻り: (per-channel [static, fsum_lhs, fsum_rhs], 対リスト (s, zD, zS, zL) の下位 K)
struct PairOut {
    stat: [f64; 3],
    fsum_lhs: [f64; 3],
    fsum_rhs: [f64; 3],
    low: Vec<(f64, [f64; 3])>,
}

fn pair_spectral(
    w1: &[f64],
    v1: &[f64],
    w2: &[f64],
    v2: &[f64],
    dim: usize,
    ops: &[Vec<f64>; 3],
    keep: usize,
) -> PairOut {
    let nocc = dim / 2;
    // M_A = V2^T O_A V1 (occ 列)
    let tv = |o: &[f64]| -> Vec<f64> {
        let mut t = vec![0.0f64; dim * nocc];
        for ccol in 0..nocc {
            for r in 0..dim {
                let mut s = 0.0;
                for k in 0..dim {
                    s += o[r + k * dim] * v1[k + ccol * dim];
                }
                t[r + ccol * dim] = s;
            }
        }
        t
    };
    let ts: Vec<Vec<f64>> = ops.iter().map(|o| tv(o)).collect();
    let mut out = PairOut {
        stat: [0.0; 3],
        fsum_lhs: [0.0; 3],
        fsum_rhs: [0.0; 3],
        low: Vec::new(),
    };
    for mu in nocc..dim {
        for nu in 0..nocc {
            let mut mm = [0.0f64; 3];
            for (a, t) in ts.iter().enumerate() {
                let mut s = 0.0;
                for k in 0..dim {
                    s += v2[k + mu * dim] * t[k + nu * dim];
                }
                mm[a] = s;
            }
            let de = w2[mu] - w1[nu];
            let s_val = de * de;
            let z = [mm[0] * mm[0], mm[1] * mm[1], mm[2] * mm[2]];
            for a in 0..3 {
                out.stat[a] += 2.0 * z[a] / de;
                out.fsum_lhs[a] += 2.0 * de * z[a];
            }
            out.low.push((s_val, z));
        }
    }
    // 下位 keep 対のみ保持 (決定的: s 昇順、タイは元の走査順)
    out.low
        .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    out.low.truncate(keep);
    // f-sum RHS: ⟨[T†,[H,T]]⟩ = tr(G1 O^T X) − tr(G2 X O^T), X = h2 O − O h1。
    // ここで h1, h2 は対角化済み表現を使う: G = V_occ V_occ^T, h = V W V^T。
    // 基底変換で O' = V2^T O V1 (dim×dim), X' = W2 O' − O' W1 (対角×行列),
    // tr(G1 O^T X) = Σ_{ν occ} (O'^T X')_{νν} = Σ_{ν occ} Σ_μ O'_{μν} X'_{μν}
    // tr(G2 X O^T) = Σ_{μ occ2} Σ_ν X'_{μν} O'_{μν}
    {
        // O' 全列 (occ + unocc) が要る — tv は occ 列のみだったので全列版
        let tv_full = |o: &[f64]| -> Vec<f64> {
            let mut t = vec![0.0f64; dim * dim];
            for ccol in 0..dim {
                for r in 0..dim {
                    let mut s = 0.0;
                    for k in 0..dim {
                        s += o[r + k * dim] * v1[k + ccol * dim];
                    }
                    t[r + ccol * dim] = s;
                }
            }
            let mut op = vec![0.0f64; dim * dim];
            for ccol in 0..dim {
                for mu in 0..dim {
                    let mut s = 0.0;
                    for k in 0..dim {
                        s += v2[k + mu * dim] * t[k + ccol * dim];
                    }
                    op[mu + ccol * dim] = s;
                }
            }
            op
        };
        for (a, o) in ops.iter().enumerate() {
            let op = tv_full(o); // O'_{μν} = ⟨μ(2)|O|ν(1)⟩ 全添字
            let mut t1 = 0.0f64; // Σ_{ν occ1} Σ_μ O'_{μν} X'_{μν}
            let mut t2 = 0.0f64; // Σ_{μ occ2} Σ_ν O'_{μν} X'_{μν}
            for nu in 0..dim {
                for mu in 0..dim {
                    let xp = (w2[mu] - w1[nu]) * op[mu + nu * dim];
                    if nu < nocc {
                        t1 += op[mu + nu * dim] * xp;
                    }
                    if mu < nocc {
                        t2 += op[mu + nu * dim] * xp;
                    }
                }
            }
            out.fsum_rhs[a] = t1 - t2;
        }
    }
    out
}

/// 走査: 固定 q = 2πj/N の全ブロック対 (決定的スレッド分割)。
struct ScanOut {
    stat: [f64; 3],
    fsum_lhs: [f64; 3],
    fsum_rhs: [f64; 3],
    low: Vec<(f64, [f64; 3])>,
    chi_ref: [f64; 3], // 独立経路: chi_pair 型 Lehmann (S0 用)
}

fn spectral_scan(n: usize, m: f64, j: usize, nthreads: usize, mutate: bool, keep: usize) -> ScanOut {
    let nb = n / 2;
    let mut rows: Vec<Option<(PairOut, [f64; 3])>> = Vec::new();
    rows.resize_with(nb, || None);
    let chunk = nb.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let jz = t * chunk + i;
                    let ckz = (2.0 * PI * jz as f64 / n as f64).cos();
                    let mut eigs: Vec<(Vec<f64>, Vec<f64>)> = Vec::with_capacity(nb);
                    for jy in 0..nb {
                        let cky = (2.0 * PI * jy as f64 / n as f64).cos();
                        let h = block_h(n, m, cky, ckz);
                        eigs.push(jacobi_eigh(&h, 4 * n));
                    }
                    let dim = 4 * n;
                    let q = 2.0 * PI * j as f64 / n as f64;
                    let mut acc = PairOut {
                        stat: [0.0; 3],
                        fsum_lhs: [0.0; 3],
                        fsum_rhs: [0.0; 3],
                        low: Vec::new(),
                    };
                    let mut chi_ref = [0.0f64; 3];
                    for jy in 0..nb {
                        let ky = 2.0 * PI * jy as f64 / n as f64;
                        let mut jt = jy + j;
                        let mut sw = false;
                        while jt >= nb {
                            jt -= nb;
                            sw = !sw;
                        }
                        let sw_eff = if mutate { false } else { sw };
                        let ops = channel_vertices(n, ky, ckz, q, sw_eff);
                        let (w1, v1) = &eigs[jy];
                        let (w2, v2) = &eigs[jt];
                        let po = pair_spectral(w1, v1, w2, v2, dim, &ops, keep);
                        for a in 0..3 {
                            acc.stat[a] += po.stat[a];
                            acc.fsum_lhs[a] += po.fsum_lhs[a];
                            acc.fsum_rhs[a] += po.fsum_rhs[a];
                        }
                        acc.low.extend(po.low);
                        // 独立経路 (chi_pair 型): tv → (mu,nu) 縮約を別実装で
                        for (a, o) in ops.iter().enumerate() {
                            let nocc = dim / 2;
                            let mut tvv = vec![0.0f64; dim * nocc];
                            for ccol in 0..nocc {
                                for r in 0..dim {
                                    let mut s = 0.0;
                                    for k in 0..dim {
                                        s += o[r + k * dim] * v1[k + ccol * dim];
                                    }
                                    tvv[r + ccol * dim] = s;
                                }
                            }
                            let mut chi = 0.0f64;
                            for mu in nocc..dim {
                                for nu in 0..nocc {
                                    let mut mr = 0.0f64;
                                    for k in 0..dim {
                                        mr += v2[k + mu * dim] * tvv[k + nu * dim];
                                    }
                                    chi += 2.0 * mr * mr / (w2[mu] - w1[nu]);
                                }
                            }
                            chi_ref[a] += chi;
                        }
                    }
                    acc.low.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                    acc.low.truncate(keep);
                    *slot = Some((acc, chi_ref));
                }
            });
        }
    });
    let vol = (n * n * n) as f64;
    let mut out = ScanOut {
        stat: [0.0; 3],
        fsum_lhs: [0.0; 3],
        fsum_rhs: [0.0; 3],
        low: Vec::new(),
        chi_ref: [0.0; 3],
    };
    for r in rows {
        let (acc, cr) = r.unwrap();
        for a in 0..3 {
            out.stat[a] += acc.stat[a] / vol;
            out.fsum_lhs[a] += acc.fsum_lhs[a] / vol;
            out.fsum_rhs[a] += acc.fsum_rhs[a] / vol;
            out.chi_ref[a] += cr[a] / vol;
        }
        out.low.extend(acc.low);
    }
    out.low.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    out.low.truncate(keep);
    out
}

// ---------------- dense 照合 (N=8) ----------------

fn build_h_dense(n: usize, m_stag: f64) -> Vec<f64> {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut h = vec![0.0f64; ns * ns];
    let add = |h: &mut Vec<f64>, i: usize, j: usize, t: f64| {
        h[j + i * ns] += t;
        h[i + j * ns] += t;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                let tw = if x == n - 1 { -1.0 } else { 1.0 };
                add(&mut h, i, idx((x + 1) % n, y, z), 0.5 * tw);
                let ey = if x % 2 == 0 { 1.0 } else { -1.0 };
                add(&mut h, i, idx(x, (y + 1) % n, z), 0.5 * ey);
                let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                add(&mut h, i, idx(x, y, (z + 1) % n), 0.5 * ez);
                let sgn = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                h[i + i * ns] += m_stag * sgn;
            }
        }
    }
    h
}

fn dense_vertex(n: usize, qy: f64, which: usize) -> (Vec<f64>, Vec<f64>) {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut re = vec![0.0f64; ns * ns];
    let mut im = vec![0.0f64; ns * ns];
    let addc = |i: usize, j2: usize, t: f64, ph: f64, re: &mut Vec<f64>, im: &mut Vec<f64>| {
        let (cp, sp) = (ph.cos(), ph.sin());
        re[j2 + i * ns] += t * cp;
        re[i + j2 * ns] += t * cp;
        im[j2 + i * ns] += t * sp;
        im[i + j2 * ns] += t * sp;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                if which == 1 {
                    let tw = if x == n - 1 { -1.0 } else { 1.0 };
                    addc(i, idx((x + 1) % n, y, z), 0.5 * tw, qy * y as f64, &mut re, &mut im);
                }
                if which == 2 {
                    let ey = if x % 2 == 0 { 1.0 } else { -1.0 };
                    addc(
                        i,
                        idx(x, (y + 1) % n, z),
                        0.5 * ey,
                        qy * (y as f64 + 0.5),
                        &mut re,
                        &mut im,
                    );
                }
                if which == 3 {
                    let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                    addc(i, idx(x, y, (z + 1) % n), 0.5 * ez, qy * y as f64, &mut re, &mut im);
                }
            }
        }
    }
    (re, im)
}

fn main() {
    self_test();
    println!(
        "=== v26.7 (II/II) v267_spectral — 有限サイズ spectral measure と pole 判定 (経路 B) ===\n"
    );
    println!("事前登録: spec §7 (157ca53)。pole 判定 = residue の体積 scaling (s_n → 0 ∧ Z_n → Z_* > 0)。");
    println!("自由場の登録予想 (v26.5): pole なし (連続体のみ)。q = 2π/16 固定, N ∈ {{16,32,64}}。\n");
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(4);
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!(
            "  [{}] {}  {}",
            if ok { "PASS" } else { "FAIL" },
            name,
            detail
        );
        if !ok {
            nfail += 1;
        }
    };
    let ms = [0.0f64, 0.5];
    let ns_list = [16usize, 32, 64];
    let q_phys = 2.0 * PI / 16.0;
    let keep = 256usize;
    let chnames = ["D (spin-2)", "S (P0s)", "L (縦=ゲージ)"];

    // ---- [S1a] dense f-sum 照合 (N=8, q = 2π/8·1) ----
    {
        let n = 8usize;
        let ns = n * n * n;
        let mut worst = 0.0f64;
        for &m in &ms {
            let h = build_h_dense(n, m);
            let (w, v) = jacobi_eigh(&h, ns);
            let qy = 2.0 * PI / n as f64;
            let r2 = (2.0f64).sqrt();
            // チャネル D の複素頂点
            let (xr, xi) = dense_vertex(n, qy, 1);
            let (zr, zi) = dense_vertex(n, qy, 3);
            let dr: Vec<f64> = xr.iter().zip(&zr).map(|(a, b)| (a - b) / r2).collect();
            let di: Vec<f64> = xi.iter().zip(&zi).map(|(a, b)| (a - b) / r2).collect();
            // LHS: 対和 Σ 2ΔE|M|²
            let nocc = ns / 2;
            let tv = |o: &[f64]| -> Vec<f64> {
                let mut t = vec![0.0f64; ns * nocc];
                for c in 0..nocc {
                    for r in 0..ns {
                        let mut s = 0.0;
                        for k in 0..ns {
                            s += o[k + r * ns] * v[k + c * ns];
                        }
                        t[r + c * ns] = s;
                    }
                }
                t
            };
            let tre = tv(&dr);
            let tim = tv(&di);
            let mut lhs = 0.0f64;
            for mu in nocc..ns {
                for nu in 0..nocc {
                    let (mut mr, mut mi) = (0.0f64, 0.0f64);
                    for k in 0..ns {
                        let vm = v[k + mu * ns];
                        mr += vm * tre[k + nu * ns];
                        mi += vm * tim[k + nu * ns];
                    }
                    lhs += 2.0 * (w[mu] - w[nu]) * (mr * mr + mi * mi);
                }
            }
            // RHS: ⟨[T†,[H,T]]⟩ — 1 体行列 [T†,[h,T]] の占有 trace (複素 T = dr + i di)
            // [h,T] = hT − Th (実部/虚部それぞれ)。[T†, X] = T†X − XT†, T† = dr^T − i di^T
            // (dr, di は対称行列として構成済み — T† は (dr − i di))
            let matmul = |a: &[f64], b: &[f64]| -> Vec<f64> {
                let mut c = vec![0.0f64; ns * ns];
                for i in 0..ns {
                    for k in 0..ns {
                        let aik = a[k + i * ns];
                        if aik == 0.0 {
                            continue;
                        }
                        for j2 in 0..ns {
                            c[j2 + i * ns] += aik * b[j2 + k * ns];
                        }
                    }
                }
                c
            };
            let sub = |a: &[f64], b: &[f64]| -> Vec<f64> {
                a.iter().zip(b).map(|(x, y)| x - y).collect()
            };
            // X = [h, T]: Xr = hTr − Trh, Xi = hTi − Tih
            let xr2 = sub(&matmul(&h, &dr), &matmul(&dr, &h));
            let xi2 = sub(&matmul(&h, &di), &matmul(&di, &h));
            // Y = [T†, X] = (Tr − iTi)(Xr + iXi) − (Xr + iXi)(Tr − iTi); 実部のみ要る
            // Re Y = TrXr + TiXi − XrTr − XiTi
            let yre = sub(
                &{
                    let mut s1 = matmul(&dr, &xr2);
                    let s2 = matmul(&di, &xi2);
                    for k in 0..ns * ns {
                        s1[k] += s2[k];
                    }
                    s1
                },
                &{
                    let mut s1 = matmul(&xr2, &dr);
                    let s2 = matmul(&xi2, &di);
                    for k in 0..ns * ns {
                        s1[k] += s2[k];
                    }
                    s1
                },
            );
            let mut rhs = 0.0f64;
            for nu in 0..nocc {
                for r in 0..ns {
                    let mut s = 0.0;
                    for k in 0..ns {
                        s += yre[k + r * ns] * v[k + nu * ns];
                    }
                    rhs += v[r + nu * ns] * s;
                }
            }
            worst = worst.max((lhs / rhs - 1.0).abs());
        }
        check(
            "[S1a] dense f-sum: Σ 2ΔE Z_D = ⟨[T_D†,[H,T_D]]⟩ (N=8, m 2 種)",
            worst < 1e-10,
            format!("max 相対 = {:.1e} ({} s)", worst, t0.elapsed().as_secs()),
        );
    }

    // ---- 走査 (q 固定の体積列) ----
    let mut scans: Vec<Vec<ScanOut>> = Vec::new(); // [mi][ni]
    for &m in &ms {
        let mut per_n = Vec::new();
        for &n in &ns_list {
            let j = n / 16;
            let s = spectral_scan(n, m, j, nthreads, false, keep);
            println!(
                "    [走査] N={} m={:.1} (j={}) 完了 ({} s) — s_min = {:.6}, χ_D = {:.6}",
                n,
                m,
                j,
                t0.elapsed().as_secs(),
                s.low[0].0,
                s.stat[0]
            );
            per_n.push(s);
        }
        scans.push(per_n);
    }

    // ---- [S0] 内部恒等 (対和 = chi_pair 型 Lehmann) ----
    {
        let mut worst = 0.0f64;
        for mi in 0..ms.len() {
            for ni in 0..ns_list.len() {
                let s = &scans[mi][ni];
                for a in 0..3 {
                    worst = worst.max((s.stat[a] - s.chi_ref[a]).abs());
                }
            }
        }
        check(
            "[S0] 内部恒等: 対和 Σ 2Z/ΔE = ブロック Lehmann (全 N, m, チャネル)",
            worst < 1e-12,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- [S1b] block f-sum (全 N) ----
    {
        let mut worst = 0.0f64;
        for mi in 0..ms.len() {
            for ni in 0..ns_list.len() {
                let s = &scans[mi][ni];
                for a in 0..3 {
                    worst = worst.max((s.fsum_lhs[a] / s.fsum_rhs[a] - 1.0).abs());
                }
            }
        }
        check(
            "[S1b] block f-sum: Σ 2ΔE Z = ⟨[T†,[H,T]]⟩ (全 N, m, チャネル)",
            worst < 1e-10,
            format!("max 相対 = {:.1e}", worst),
        );
    }

    // ---- [S2] 正値性 ----
    {
        let mut min_z = f64::INFINITY;
        let mut min_s = f64::INFINITY;
        for mi in 0..ms.len() {
            for ni in 0..ns_list.len() {
                for &(s, z) in &scans[mi][ni].low {
                    min_s = min_s.min(s);
                    for a in 0..3 {
                        min_z = min_z.min(z[a]);
                    }
                }
            }
        }
        check(
            "[S2] 正値性: Z ≥ 0 (構成的) かつ s_min > 0",
            min_z >= 0.0 && min_s > 0.0,
            format!("min Z = {:.1e}, min s = {:.4}", min_z, min_s),
        );
    }

    // ---- [S3] threshold (m=0.5): s ≥ 4m² ----
    {
        let mut worst = f64::INFINITY;
        for ni in 0..ns_list.len() {
            worst = worst.min(scans[1][ni].low[0].0);
        }
        check(
            "[S3] threshold: s_min ≥ 4m² = 1.0 (m=0.5 — E(k) ≥ m の厳密帰結)",
            worst >= 1.0 - 1e-12,
            format!("min s = {:.6} ≥ 1.0", worst),
        );
    }

    // ---- [S4] massless edge (記録) ----
    {
        println!("    [S4 表] m=0 の端 s_min/q²: ");
        for ni in 0..ns_list.len() {
            let s = &scans[0][ni];
            println!(
                "      N={}: s_min = {:.6} (s_min/q² = {:.4})",
                ns_list[ni],
                s.low[0].0,
                s.low[0].0 / (q_phys * q_phys)
            );
        }
        println!("      — 連続体の光円錐端は s = q² = {:.6} (branch 記録)", q_phys * q_phys);
    }

    // ---- [S5] pole 判定 (branch — 主結果) ----
    {
        // 最低クラスタ (s_min の相対 1e-9 以内) の Σ Z / V
        let w1 = |mi: usize, ni: usize, a: usize| -> f64 {
            let s = &scans[mi][ni];
            let smin = s.low[0].0;
            let vol = (ns_list[ni].pow(3)) as f64;
            s.low
                .iter()
                .filter(|&&(sv, _)| sv < smin * (1.0 + 1e-9))
                .map(|&(_, z)| z[a])
                .sum::<f64>()
                / vol
        };
        println!("    [S5 表] W₁(N) = 最低クラスタ Σ Z / V (チャネル D / L):");
        let mut resolved = true;
        for (mi, &m) in ms.iter().enumerate() {
            let wd: Vec<f64> = (0..3).map(|ni| w1(mi, ni, 0)).collect();
            let wl: Vec<f64> = (0..3).map(|ni| w1(mi, ni, 2)).collect();
            println!(
                "      m={:.1}: D: {:.3e} → {:.3e} → {:.3e} / L: {:.3e} → {:.3e} → {:.3e}",
                m, wd[0], wd[1], wd[2], wl[0], wl[1], wl[2]
            );
            let mono_dec = wd[0] > wd[1] && wd[1] > wd[2];
            let mono_inc = wd[0] < wd[1] && wd[1] < wd[2];
            if !(mono_dec || mono_inc) {
                resolved = false;
            }
            if mono_dec {
                println!(
                    "        ⇒ branch α (m={:.1}): W₁ は体積とともに消える — **pole なし** (残差冪 ~ N^{:.2})",
                    m,
                    (wd[2] / wd[0]).ln() / (ns_list[2] as f64 / ns_list[0] as f64).ln()
                );
            } else if mono_inc {
                println!("        ⇒ branch β (m={:.1}): W₁ が成長 — pole 様 (登録予想と不一致!)", m);
            }
        }
        check(
            "[S5] pole 判定の分解能: W₁(N) の順序が単調に確定 (branch α/β の判定可能性)",
            resolved,
            "branch は上の欄が一次ソース".to_string(),
        );
    }

    // ---- [S6] 変異検出 ----
    {
        let n = 16usize;
        let good = &scans[0][0];
        let bad = spectral_scan(n, 0.0, 1, nthreads, true, keep);
        let dev = (bad.stat[0] - good.stat[0])
            .abs()
            .max((bad.fsum_lhs[0] - good.fsum_lhs[0]).abs());
        check(
            "[S6] 変異: 折返しスワップ落とし → S0/S1 の量が変化 (> 1e-4)",
            dev > 1e-4,
            format!("逸脱 {:.2e}", dev),
        );
    }

    // ---- スペクトル表 (定量の一次ソース) ----
    println!("\n    [spectral 表 (N=64)] チャネル | m | χ (=Σ2Z/ΔE) | f-sum (=Σ2ΔE·Z) | s_min | W₁");
    for (mi, &m) in ms.iter().enumerate() {
        for a in 0..3 {
            let s = &scans[mi][2];
            let vol = (ns_list[2].pow(3)) as f64;
            let smin = s.low[0].0;
            let w1: f64 = s
                .low
                .iter()
                .filter(|&&(sv, _)| sv < smin * (1.0 + 1e-9))
                .map(|&(_, z)| z[a])
                .sum::<f64>()
                / vol;
            println!(
                "      {} m={:.1}: {:.6} | {:.6} | {:.6} | {:.3e}",
                chnames[a], m, s.stat[a], s.fsum_lhs[a], smin, w1
            );
        }
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.7-II".into())),
        ("kind".into(), Json::Str("spectral_measure_pole_search".into())),
        ("q_phys".into(), Json::Num(q_phys)),
        (
            "w1_scaling".into(),
            Json::Arr(
                (0..2)
                    .map(|mi| {
                        let w: Vec<Json> = (0..3)
                            .map(|ni| {
                                let s = &scans[mi][ni];
                                let smin = s.low[0].0;
                                let vol = (ns_list[ni].pow(3)) as f64;
                                Json::Num(
                                    s.low
                                        .iter()
                                        .filter(|&&(sv, _)| sv < smin * (1.0 + 1e-9))
                                        .map(|&(_, z)| z[0])
                                        .sum::<f64>()
                                        / vol,
                                )
                            })
                            .collect();
                        Json::Obj(vec![
                            ("m".into(), Json::Num(ms[mi])),
                            ("w1_D_by_n".into(), Json::Arr(w)),
                            (
                                "s_min_by_n".into(),
                                Json::Arr((0..3).map(|ni| Json::Num(scans[mi][ni].low[0].0)).collect()),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v267_spectral.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **spectral 器械認証が成立、S5 の branch が主結果** — 解釈は docs/uft-v26.7.md §II へ"
        } else {
            "FAIL — 分岐 (b) 対分解の誤り / (c) 体積列の拡張"
        }
    );
    println!(
        "\n総合判定: {} ({} s)",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" },
        t0.elapsed().as_secs()
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
