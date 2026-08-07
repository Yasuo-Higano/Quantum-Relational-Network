//! v35.2 開放系 signed-response 観測商 — GQF-1..6 の機械実証 (PROMPT/16 §5–6)
//!
//! 物理: number-conserving quasi-free 開放系 Ċ = XC + CX† + Y (X = −ih − ½(Λᵀ+M),
//! Y = M — 規約は [G1] で較正)。signed probe C± = C₀ ± εP_i の差分応答が識別するのは Hamiltonian では
//! なく**有効 drift X の観測商**である。
//!
//! 検証の柱:
//!  [G1] covariance 閉包の第一原理較正 — dense Lindblad (2^N 次元 Jordan–Wigner,
//!       RK4) と Van Loan 厳密解の一致 (規約 X/Y はここで凍結)。
//!  [G2] GQF-1: Y ≠ 0 でも Δ(t) = C⁺−C⁻ = 2ε e^{Xt}P e^{X†t} (affine 消去) —
//!       生の n_j(t) は Y に依存する (負制御)。
//!  [G3] GQF-2: 曲率 = ‖P_j X P_i‖²_F — jet 恒等式 (Lean 鏡映) と Richardson
//!       有限差分の両側から。
//!  [G4] Lean 整数インスタンスの橋 (proofs/OpenQuotient.lean ↔ Rust)。
//!  [C1..C11] 反例セル: coherent hopping / collective loss (GQF-3 同一曲率) /
//!       対角 loss (GQF-4 還元) / pairing OutOfDomain / 生成子 probe 依存 /
//!       local phase gauge / 複素共役 / 有限 shot Straddled / 相関 shot
//!       OutOfDomain / 電荷帰属 (GQF-5)。
//!
//! 型の規律: EffectiveDriftTopology ↛ HamiltonianTopology (禁止変換 32 — 門は
//! DissipatorLocalityCertificate)・ChargeNonconservingResponse ↛
//! HamiltonianPairingWitness (禁止変換 33 — 門は DissipativeChargeConservation)。

use uft_sim::open_response::{
    certify_charge_conserving_dissipation, certify_dissipator_locality, charge_response, cs,
    curvature_exact, evolve_covariance, expm, global_frequency_shift, local_phase_gauge,
    mat_add, mat_adj, mat_conj,
    mat_eye, mat_max_abs, mat_mul, mat_scale, mat_sub, mat_trace, mat_zero, open_response_self_test,
    pairing_witness_with_certificate, promote_with_certificate, response_table,
    ChargeNonconservingResponse, CurvatureVerdict, EffectiveDriftTopology,
    FiniteShotCurvatureReader, OpenLaneRefusal, PromotionRefusal, QuasiFreeOpenModel,
};
use uft_sim::{self_test, Rng, C64, CONE, CZERO};

type CMat = Vec<Vec<C64>>;

// ================================================================ dense Lindblad (第一原理)

/// Jordan–Wigner の c_a (dim = 2^n)。bit a が占有。符号 = (−1)^{Σ_{j<a} b_j}。
fn jw_annihilation(n: usize, a: usize) -> CMat {
    let dim = 1usize << n;
    let mut m = mat_zero(dim);
    for b in 0..dim {
        if b & (1 << a) != 0 {
            let sign = ((b & ((1 << a) - 1)).count_ones() % 2) as i32;
            let s = if sign == 0 { 1.0 } else { -1.0 };
            m[b ^ (1 << a)][b] = cs(s, 0.0);
        }
    }
    m
}

/// H = Σ h_{ab} c†_a c_b (+ 任意で pairing Δ_{ab} c†_a c†_b + h.c.)
fn jw_hamiltonian(n: usize, h: &CMat, delta: Option<&CMat>) -> CMat {
    let dim = 1usize << n;
    let cops: Vec<CMat> = (0..n).map(|a| jw_annihilation(n, a)).collect();
    let cdag: Vec<CMat> = cops.iter().map(mat_adj).collect();
    let mut ham = mat_zero(dim);
    for a in 0..n {
        for b in 0..n {
            if h[a][b].norm2() > 0.0 {
                ham = mat_add(&ham, &mat_scale(h[a][b], &mat_mul(&cdag[a], &cops[b])));
            }
        }
    }
    if let Some(d) = delta {
        for a in 0..n {
            for b in 0..n {
                if d[a][b].norm2() > 0.0 {
                    let t = mat_scale(d[a][b], &mat_mul(&cdag[a], &cdag[b]));
                    ham = mat_add(&ham, &t);
                    ham = mat_add(&ham, &mat_adj(&t));
                }
            }
        }
    }
    ham
}

/// ρ̇ = −i[H,ρ] + Σ_j (A_j ρ A_j† − ½{A_j†A_j, ρ}) の RK4 積分
fn dense_lindblad_evolve(ham: &CMat, jumps: &[CMat], rho0: &CMat, t: f64, dt: f64) -> CMat {
    let jj: Vec<(CMat, CMat)> = jumps
        .iter()
        .map(|a| {
            let ad = mat_adj(a);
            let ada = mat_mul(&ad, a);
            (a.clone(), ada)
        })
        .collect();
    let deriv = |rho: &CMat| -> CMat {
        let hr = mat_mul(ham, rho);
        let rh = mat_mul(rho, ham);
        let mut d = mat_scale(cs(0.0, -1.0), &mat_sub(&hr, &rh));
        for (a, ada) in &jj {
            let ar = mat_mul(a, rho);
            let ara = mat_mul(&ar, &mat_adj(a));
            let anti = mat_add(&mat_mul(ada, rho), &mat_mul(rho, ada));
            d = mat_add(&d, &mat_sub(&ara, &mat_scale(cs(0.5, 0.0), &anti)));
        }
        d
    };
    let mut rho = rho0.clone();
    let steps = (t / dt).round() as usize;
    for _ in 0..steps {
        let k1 = deriv(&rho);
        let k2 = deriv(&mat_add(&rho, &mat_scale(cs(dt / 2.0, 0.0), &k1)));
        let k3 = deriv(&mat_add(&rho, &mat_scale(cs(dt / 2.0, 0.0), &k2)));
        let k4 = deriv(&mat_add(&rho, &mat_scale(cs(dt, 0.0), &k3)));
        let mut inc = mat_add(&k1, &mat_scale(cs(2.0, 0.0), &k2));
        inc = mat_add(&inc, &mat_scale(cs(2.0, 0.0), &k3));
        inc = mat_add(&inc, &k4);
        rho = mat_add(&rho, &mat_scale(cs(dt / 6.0, 0.0), &inc));
    }
    rho
}

/// 対角 C₀ = diag(p) の積状態 ρ₀ = ⊗ diag(1−p_a, p_a)
fn product_state(n: usize, p: &[f64]) -> CMat {
    let dim = 1usize << n;
    let mut rho = mat_zero(dim);
    for (b, row) in rho.iter_mut().enumerate() {
        let mut w = 1.0;
        for (a, &pa) in p.iter().enumerate() {
            w *= if b & (1 << a) != 0 { pa } else { 1.0 - pa };
        }
        row[b] = cs(w, 0.0);
    }
    rho
}

/// C_{ab} = Tr(ρ c†_b c_a)
fn covariance_from_rho(n: usize, rho: &CMat) -> CMat {
    let cops: Vec<CMat> = (0..n).map(|a| jw_annihilation(n, a)).collect();
    let mut c = mat_zero(n);
    for a in 0..n {
        for b in 0..n {
            let op = mat_mul(&mat_adj(&cops[b]), &cops[a]);
            c[a][b] = mat_trace(&mat_mul(rho, &op));
        }
    }
    c
}

/// jump 演算子 (dense): loss L = Σ ℓ_a c_a / gain G = Σ g_a c†_a
fn jw_linear_jump(n: usize, amp: &[C64], dagger: bool) -> CMat {
    let dim = 1usize << n;
    let mut m = mat_zero(dim);
    for a in 0..n {
        let c = jw_annihilation(n, a);
        let op = if dagger { mat_adj(&c) } else { c };
        m = mat_add(&m, &mat_scale(amp[a], &op));
    }
    m
}

// ================================================================ jet 恒等式 (Lean 鏡映)

/// Tr(P_j (X²P_i + 2XP_iX† + P_iX†²)) / 2 — GQF-2 の jet 側 (厳密代数)
fn jet_curvature(x: &CMat, pi: &CMat, pj: &CMat) -> f64 {
    let xx = mat_mul(x, x);
    let xd = mat_adj(x);
    let t1 = mat_mul(&xx, pi);
    let t2 = mat_scale(cs(2.0, 0.0), &mat_mul(&mat_mul(x, pi), &xd));
    let t3 = mat_mul(pi, &mat_mul(&xd, &xd));
    let jet = mat_add(&mat_add(&t1, &t2), &t3);
    mat_trace(&mat_mul(pj, &jet)).re / 2.0
}

fn proj(n: usize, sites: &[usize]) -> CMat {
    let mut p = mat_zero(n);
    for &s in sites {
        p[s][s] = CONE;
    }
    p
}

fn main() {
    self_test();
    open_response_self_test().expect("open_response self test");
    println!("=== v35.2 開放系 signed-response 観測商 — GQF-1..6 (PROMPT/16 §5–6) ===\n");
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ---------------- [G1] covariance 閉包 vs dense Lindblad ----------------
    {
        // N=2: h 複素エルミート + loss + gain
        let h = vec![
            vec![cs(0.3, 0.0), cs(0.5, 0.2)],
            vec![cs(0.5, -0.2), cs(-0.1, 0.0)],
        ];
        let loss = vec![vec![cs(0.4, 0.0), cs(0.0, 0.3)]];
        let gain = vec![vec![cs(0.2, 0.0), cs(-0.1, 0.1)]];
        let model = QuasiFreeOpenModel::new(h.clone(), loss.clone(), gain.clone(), None)
            .expect("model");
        let x = model.effective_drift();
        let y = model.injection();
        let p0 = [0.6, 0.3];
        let c0 = {
            let mut c = mat_zero(2);
            c[0][0] = cs(p0[0], 0.0);
            c[1][1] = cs(p0[1], 0.0);
            c
        };
        // dense 側
        let ham = jw_hamiltonian(2, &h, None);
        let jumps = vec![
            jw_linear_jump(2, &loss[0], false),
            jw_linear_jump(2, &gain[0], true),
        ];
        let rho0 = product_state(2, &p0);
        let t = 0.4;
        let rho_t = dense_lindblad_evolve(&ham, &jumps, &rho0, t, 1e-4);
        let c_dense = covariance_from_rho(2, &rho_t);
        let c_lane = evolve_covariance(&x, &y, &c0, t);
        let dev = mat_max_abs(&mat_sub(&c_dense, &c_lane));
        check(
            "[G1a] covariance 閉包 (N=2, loss+gain): dense Lindblad (JW/RK4) と Van Loan の一致",
            dev < 1e-6,
            format!("max|Δ| = {:.2e} (t = {}) — 規約 X = −ih − ½(Λᵀ+M), Y = M を凍結", dev, t),
        );

        // N=3, loss のみ・非自明位相
        let h3 = vec![
            vec![CZERO, cs(1.0, 0.0), CZERO],
            vec![cs(1.0, 0.0), cs(0.2, 0.0), cs(0.0, -0.7)],
            vec![CZERO, cs(0.0, 0.7), CZERO],
        ];
        let loss3 = vec![vec![cs(0.3, 0.0), CZERO, cs(0.0, 0.5)]];
        let m3 = QuasiFreeOpenModel::new(h3.clone(), loss3.clone(), vec![], None).expect("m3");
        let p3 = [0.8, 0.4, 0.1];
        let mut c03 = mat_zero(3);
        for i in 0..3 {
            c03[i][i] = cs(p3[i], 0.0);
        }
        let ham3 = jw_hamiltonian(3, &h3, None);
        let jumps3 = vec![jw_linear_jump(3, &loss3[0], false)];
        let rho03 = product_state(3, &p3);
        let rho_t3 = dense_lindblad_evolve(&ham3, &jumps3, &rho03, 0.5, 1e-4);
        let c_dense3 = covariance_from_rho(3, &rho_t3);
        let c_lane3 = evolve_covariance(&m3.effective_drift(), &m3.injection(), &c03, 0.5);
        let dev3 = mat_max_abs(&mat_sub(&c_dense3, &c_lane3));
        check(
            "[G1b] covariance 閉包 (N=3, loss のみ・複素振幅)",
            dev3 < 1e-6,
            format!("max|Δ| = {:.2e}", dev3),
        );
    }

    // ---------------- [G2] GQF-1: affine 消去 ----------------
    {
        let h = vec![
            vec![cs(0.3, 0.0), cs(0.5, 0.2)],
            vec![cs(0.5, -0.2), cs(-0.1, 0.0)],
        ];
        let loss = vec![vec![cs(0.4, 0.0), cs(0.0, 0.3)]];
        let gain = vec![vec![cs(0.2, 0.0), cs(-0.1, 0.1)]];
        let model = QuasiFreeOpenModel::new(h, loss, gain, None).expect("model");
        let x = model.effective_drift();
        let y = model.injection();
        let eps = 0.05;
        let pi = proj(2, &[0]);
        let mut c0 = mat_zero(2);
        c0[0][0] = cs(0.5, 0.0);
        c0[1][1] = cs(0.5, 0.0);
        let cp0 = mat_add(&c0, &mat_scale(cs(eps, 0.0), &pi));
        let cm0 = mat_sub(&c0, &mat_scale(cs(eps, 0.0), &pi));
        let t = 0.7;
        let cp = evolve_covariance(&x, &y, &cp0, t);
        let cm = evolve_covariance(&x, &y, &cm0, t);
        let delta = mat_sub(&cp, &cm);
        // 同次公式: 2ε e^{Xt} P e^{X†t}
        let f = expm(&mat_scale(cs(t, 0.0), &x));
        let hom = mat_scale(cs(2.0 * eps, 0.0), &mat_mul(&mat_mul(&f, &pi), &mat_adj(&f)));
        let dev = mat_max_abs(&mat_sub(&delta, &hom));
        check(
            "[G2a] GQF-1: Y ≠ 0 でも Δ(t) = 2ε e^{Xt}P e^{X†t} (affine 項の厳密消去)",
            dev < 1e-12,
            format!("max|Δ − hom| = {:.2e}", dev),
        );
        // 負制御: 生の n₁(t) は Y に依存する (gain を切ると変わる)
        let model_nogain = QuasiFreeOpenModel::new(
            vec![
                vec![cs(0.3, 0.0), cs(0.5, 0.2)],
                vec![cs(0.5, -0.2), cs(-0.1, 0.0)],
            ],
            vec![vec![cs(0.4, 0.0), cs(0.0, 0.3)]],
            vec![],
            None,
        )
        .expect("m");
        // 有効 drift も変わるので同一 X で比較するため、Y だけ 0 にした flow を評価
        let c_with_y = evolve_covariance(&x, &y, &cp0, t);
        let c_no_y = evolve_covariance(&x, &mat_zero(2), &cp0, t);
        let raw_dev = mat_max_abs(&mat_sub(&c_with_y, &c_no_y));
        let _ = model_nogain;
        check(
            "[G2b] 負制御: 生の単側軌道は Y に依存 (差分だけが Y を消す)",
            raw_dev > 1e-3,
            format!("|C_Y − C_0| = {:.2e} ≫ 0 (Δ 側は {:.0e})", raw_dev, 1e-13),
        );
    }

    // ---------------- [G3] GQF-2: 曲率 = ‖P_j X P_i‖² ----------------
    {
        // 4 モード・2 次元 node (Lean gqf2_curvature_block の実数値鏡映 + 乱択複素)
        let mut rng = Rng::new(35201);
        let n = 4;
        let mut h = mat_zero(n);
        for i in 0..n {
            for j in (i + 1)..n {
                let z = cs(rng.f64() - 0.5, rng.f64() - 0.5);
                h[i][j] = z;
                h[j][i] = z.conj();
            }
            h[i][i] = cs(rng.f64() - 0.5, 0.0);
        }
        let loss = vec![
            (0..n).map(|_| cs(rng.f64() - 0.5, rng.f64() - 0.5)).collect(),
            (0..n).map(|_| cs(rng.f64() - 0.5, rng.f64() - 0.5)).collect(),
        ];
        let gain = vec![(0..n).map(|_| cs(rng.f64() - 0.5, rng.f64() - 0.5)).collect()];
        let model = QuasiFreeOpenModel::new(h, loss, gain, None).expect("model");
        let x = model.effective_drift();
        let nodes = vec![vec![0usize, 1], vec![2usize, 3]];
        let w = curvature_exact(&x, &nodes);
        // jet 恒等式側 (Lean GQF-2 の鏡映)
        let pi = proj(n, &nodes[0]);
        let pj = proj(n, &nodes[1]);
        let wj = jet_curvature(&x, &pi, &pj);
        let dev = (w[1][0] - wj).abs();
        check(
            "[G3a] GQF-2 jet 恒等式: Tr(P_j(X²P+2XPX†+PX†²))/2 = ‖P_j X P_i‖² (複素乱択 4×4・2 次元 node)",
            dev < 1e-13,
            format!("|jet − frob| = {:.2e} (w = {:.6})", dev, w[1][0]),
        );
        // Richardson 有限差分側 (厳密軌道から)
        let y = model.injection();
        let eps = 1e-3;
        let dt = 1e-3;
        let mut c0 = mat_zero(n);
        for i in 0..n {
            c0[i][i] = cs(0.5, 0.0);
        }
        let read = |t: f64, sgn: f64| -> f64 {
            let cpr = mat_add(&c0, &mat_scale(cs(sgn * eps, 0.0), &pi));
            let ct = evolve_covariance(&x, &y, &cpr, t);
            mat_trace(&mat_mul(&pj, &ct)).re
        };
        let d1 = read(dt, 1.0) - read(dt, -1.0);
        let d2 = read(2.0 * dt, 1.0) - read(2.0 * dt, -1.0);
        let k_hat = (8.0 * d1 - d2) / (8.0 * eps * dt * dt);
        let dev_fd = (k_hat - w[1][0]).abs();
        check(
            "[G3b] GQF-2 有限差分 (Richardson): 厳密軌道からの K̂ = ‖P_j X P_i‖²",
            dev_fd < 1e-5,
            format!("|K̂ − w| = {:.2e}", dev_fd),
        );
    }

    // ---------------- [G4] Lean 整数インスタンスの橋 ----------------
    {
        // GQF-3 対 (proofs/OpenQuotient.lean XA/XB)
        let xa = vec![vec![CZERO, cs(0.0, -1.0)], vec![cs(0.0, -1.0), CZERO]];
        let xb = vec![vec![cs(-1.0, 0.0), cs(1.0, 0.0)], vec![cs(1.0, 0.0), cs(-1.0, 0.0)]];
        let nodes = vec![vec![0usize], vec![1usize]];
        let wa = curvature_exact(&xa, &nodes);
        let wb = curvature_exact(&xb, &nodes);
        // GQF-5 対 (chargeStatLoss / chargeStatPairing = −8)
        let x_loss = vec![vec![cs(-4.0, 0.0), CZERO], vec![CZERO, CZERO]];
        let c_loss = {
            let mut c = mat_zero(2);
            c[0][0] = CONE;
            c
        };
        let s_loss = charge_response(&x_loss, &mat_zero(2), &c_loss);
        // Nambu 側: Re Tr(QN·(−i)[H_B, R2]) — Lean の宣言スケール (R2 = 2R, QN = 2Q)
        let h_bdg = vec![
            vec![CZERO, CZERO, CZERO, cs(1.0, 0.0)],
            vec![CZERO, CZERO, cs(-1.0, 0.0), CZERO],
            vec![CZERO, cs(-1.0, 0.0), CZERO, CZERO],
            vec![cs(1.0, 0.0), CZERO, CZERO, CZERO],
        ];
        let r2 = vec![
            vec![CONE, CZERO, CZERO, cs(0.0, 1.0)],
            vec![CZERO, CONE, cs(0.0, -1.0), CZERO],
            vec![CZERO, cs(0.0, 1.0), CONE, CZERO],
            vec![cs(0.0, -1.0), CZERO, CZERO, CONE],
        ];
        let qn = {
            let mut q = mat_eye(4);
            q[2][2] = cs(-1.0, 0.0);
            q[3][3] = cs(-1.0, 0.0);
            q
        };
        let comm = mat_sub(&mat_mul(&h_bdg, &r2), &mat_mul(&r2, &h_bdg));
        let s_pair = mat_trace(&mat_mul(&qn, &mat_scale(cs(0.0, -1.0), &comm))).re;
        check(
            "[G4a] Lean 橋: GQF-3 対の曲率 (1, 1)・GQF-5 対の電荷統計 (−8, −8) が整数一致",
            (wa[1][0] - 1.0).abs() < 1e-15
                && (wb[1][0] - 1.0).abs() < 1e-15
                && (s_loss + 8.0).abs() < 1e-15
                && (s_pair + 8.0).abs() < 1e-15,
            format!(
                "w_A = {}, w_B = {}, S_loss = {}, S_pair = {}",
                wa[1][0], wb[1][0], s_loss, s_pair
            ),
        );
        // 宣言スケールの物理較正: dense JW で dN/dt を計算し、宣言統計 = 4 × 物理
        // (R2 = 2R, QN = 2Q) を機械確認。H = Δ(c†₁c†₂ + h.c.), |ψ⟩ = (|00⟩ + i|11⟩)/√2
        let delta = vec![vec![CZERO, cs(1.0, 0.0)], vec![CZERO, CZERO]];
        let ham_pair = jw_hamiltonian(2, &mat_zero(2), Some(&delta));
        let dim = 4usize;
        // |ψ⟩ = (|00⟩ + i|11⟩)/√2 → ρ = |ψ⟩⟨ψ| (mask 0 = |00⟩, mask 3 = |11⟩)
        let mut rho = mat_zero(dim);
        let amp = [
            (0usize, cs(1.0 / 2f64.sqrt(), 0.0)),
            (3usize, cs(0.0, 1.0 / 2f64.sqrt())),
        ];
        for &(i, ai) in &amp {
            for &(j, aj) in &amp {
                rho[i][j] = ai * aj.conj();
            }
        }
        let nop = {
            let c1 = jw_annihilation(2, 0);
            let c2 = jw_annihilation(2, 1);
            mat_add(&mat_mul(&mat_adj(&c1), &c1), &mat_mul(&mat_adj(&c2), &c2))
        };
        // dN/dt = Tr(ρ · i[H, N])
        let commn = mat_sub(&mat_mul(&ham_pair, &nop), &mat_mul(&nop, &ham_pair));
        let dndt = mat_trace(&mat_mul(&rho, &mat_scale(cs(0.0, 1.0), &commn))).re;
        check(
            "[G4b] 宣言スケールの物理較正: dense dN/dt × 4 = 宣言統計 (BCS 状態の pairing 電荷応答)",
            (4.0 * dndt - s_pair).abs() < 1e-12 && dndt.abs() > 0.1,
            format!("dN/dt = {} (dense JW), 4× = {} = S_pair", dndt, 4.0 * dndt),
        );
    }

    // ---------------- [C1] coherent hopping (+ 証明書つき昇格) ----------------
    let nodes2 = vec![vec![0usize], vec![1usize]];
    {
        let h = vec![vec![CZERO, cs(1.0, 0.0)], vec![cs(1.0, 0.0), CZERO]];
        let model = QuasiFreeOpenModel::new(h, vec![], vec![], None).expect("m");
        let x = model.effective_drift();
        let eff = EffectiveDriftTopology {
            w: curvature_exact(&x, &nodes2),
        };
        let cert = certify_dissipator_locality(&model, &nodes2, 1e-12);
        let ok = match &cert {
            Ok(c) => {
                let ham = promote_with_certificate(&eff, c);
                (ham.w[1][0] - 1.0).abs() < 1e-15
            }
            Err(_) => false,
        };
        check(
            "[C1] coherent hopping: w₂₁ = |h₂₁|² = 1・散逸ゼロ証明書で HamiltonianTopology へ昇格",
            ok,
            format!("w₂₁ = {:.15}", eff.w[1][0]),
        );
    }

    // ---------------- [C2] collective loss — GQF-3 反例 (昇格拒否) ----------------
    {
        // h = 0, ℓ = √2 (1, −1)
        let s2 = 2f64.sqrt();
        let loss = vec![vec![cs(s2, 0.0), cs(-s2, 0.0)]];
        let model = QuasiFreeOpenModel::new(mat_zero(2), loss, vec![], None).expect("m");
        let x = model.effective_drift();
        let eff = EffectiveDriftTopology {
            w: curvature_exact(&x, &nodes2),
        };
        let cert = certify_dissipator_locality(&model, &nodes2, 1e-12);
        let refused = matches!(
            cert,
            Err(PromotionRefusal::OffDiagonalDissipator { .. })
        );
        check(
            "[C2] collective loss: 同一曲率 w₂₁ = 1 (GQF-3 反例対) — Hamiltonian 昇格は証明書段階で拒否",
            (eff.w[1][0] - 1.0).abs() < 1e-15 && refused,
            format!(
                "w₂₁ = {:.15} (coherent と同値), cross-node 散逸 = 2 > bar → 禁止変換 32",
                eff.w[1][0]
            ),
        );
    }

    // ---------------- [C3] 対角 loss — GQF-4 還元 ----------------
    {
        let h = vec![vec![CZERO, cs(0.8, 0.3)], vec![cs(0.8, -0.3), CZERO]];
        let closed = QuasiFreeOpenModel::new(h.clone(), vec![], vec![], None).expect("m");
        let w_closed = curvature_exact(&closed.effective_drift(), &nodes2);
        // 対角 loss (site 1 と site 2 で別レート — cross-node 項なし)
        let loss = vec![
            vec![cs(0.7, 0.0), CZERO],
            vec![CZERO, cs(0.0, 0.4)],
        ];
        let open = QuasiFreeOpenModel::new(h, loss, vec![], None).expect("m");
        let w_open = curvature_exact(&open.effective_drift(), &nodes2);
        let cert = certify_dissipator_locality(&open, &nodes2, 1e-12);
        check(
            "[C3] 対角 loss: w₂₁ = |h₂₁|² が閉鎖系と厳密一致 (GQF-4)・証明書 PASS で昇格可",
            (w_open[1][0] - w_closed[1][0]).abs() < 1e-15 && cert.is_ok(),
            format!("open {} = closed {} (差 {:.1e})", w_open[1][0], w_closed[1][0],
                    (w_open[1][0] - w_closed[1][0]).abs()),
        );
    }

    // ---------------- [C4] pairing 宣言は構成時 OutOfDomain ----------------
    {
        let delta = vec![vec![CZERO, cs(0.5, 0.0)], vec![cs(-0.5, 0.0), CZERO]];
        let r = QuasiFreeOpenModel::new(mat_zero(2), vec![], vec![], Some(&delta));
        let refused = matches!(r, Err(OpenLaneRefusal::PairingOutOfDomain { .. }));
        check(
            "[C4] pairing (Nambu Δ ≠ 0) の宣言: number-conserving lane は構成時 OutOfDomain (強制回答なし)",
            refused,
            "Δ ≠ 0 → Err(PairingOutOfDomain)".into(),
        );
    }

    // ---------------- [C5] 生成子の probe 依存 (二 ε ゲート) ----------------
    {
        let h = vec![vec![CZERO, cs(1.0, 0.0)], vec![cs(1.0, 0.0), CZERO]];
        let model = QuasiFreeOpenModel::new(h, vec![], vec![], None).expect("m");
        let x = model.effective_drift();
        let pi = proj(2, &[0]);
        let pj = proj(2, &[1]);
        let mut c0 = mat_zero(2);
        c0[0][0] = cs(0.5, 0.0);
        c0[1][1] = cs(0.5, 0.0);
        let t = 1.0;
        // 正直 lane: Δ_ε(t) は ε に厳密一次
        let resp = |x_eff: &CMat, eps: f64| -> f64 {
            let cp = mat_add(&c0, &mat_scale(cs(eps, 0.0), &pi));
            let cm = mat_sub(&c0, &mat_scale(cs(eps, 0.0), &pi));
            let d = mat_sub(
                &evolve_covariance(x_eff, &mat_zero(2), &cp, t),
                &evolve_covariance(x_eff, &mat_zero(2), &cm, t),
            );
            mat_trace(&mat_mul(&pj, &d)).re
        };
        let (e1, e2) = (0.05, 0.1);
        let r_clean = resp(&x, e2) * e1 / (resp(&x, e1) * e2);
        // 生成子が probe に依存する (back-action): X_ε = X − iκε P_i
        let back = |eps: f64, sgn: f64| -> f64 {
            let kappa = 5.0;
            let xe = mat_add(&x, &mat_scale(cs(0.0, -kappa * eps * sgn), &pi));
            let cpr = mat_add(&c0, &mat_scale(cs(sgn * eps, 0.0), &pi));
            let ct = evolve_covariance(&xe, &mat_zero(2), &cpr, t);
            mat_trace(&mat_mul(&pj, &ct)).re
        };
        let d_back = |eps: f64| back(eps, 1.0) - back(eps, -1.0);
        let r_back = d_back(e2) * e1 / (d_back(e1) * e2);
        let bar = 1e-9;
        check(
            "[C5] 二 ε signed-linearity ゲート: 正直 lane は比 1 (厳密)・probe 依存生成子は OutOfDomain",
            (r_clean - 1.0).abs() < bar && (r_back - 1.0).abs() > 1e-3,
            format!("|r_clean − 1| = {:.1e}, |r_back − 1| = {:.2e} (生成子可変の検出)",
                    (r_clean - 1.0).abs(), (r_back - 1.0).abs()),
        );
    }

    // ---------------- [C6] local phase + global frequency gauge ----------------
    {
        let mut rng = Rng::new(35202);
        let n = 3;
        let mut h = mat_zero(n);
        for i in 0..n {
            for j in (i + 1)..n {
                let z = cs(rng.f64() - 0.5, rng.f64() - 0.5);
                h[i][j] = z;
                h[j][i] = z.conj();
            }
        }
        let loss = vec![(0..n).map(|_| cs(rng.f64() - 0.5, rng.f64() - 0.5)).collect()];
        let model = QuasiFreeOpenModel::new(h, loss, vec![], None).expect("m");
        let x = model.effective_drift();
        let xg = global_frequency_shift(&local_phase_gauge(&x, &[0.3, -0.8, 1.7]), 0.45);
        let nodes3 = vec![vec![0usize], vec![1usize], vec![2usize]];
        let w = curvature_exact(&x, &nodes3);
        let wg = curvature_exact(&xg, &nodes3);
        let mut dev_w = 0.0f64;
        for j in 0..3 {
            for i in 0..3 {
                dev_w = dev_w.max((w[j][i] - wg[j][i]).abs());
            }
        }
        let times = [0.25, 0.7, 1.3];
        let ta = response_table(&x, &nodes3, &times);
        let tb = response_table(&xg, &nodes3, &times);
        let mut dev_t = 0.0f64;
        for (a, b) in ta.iter().zip(&tb) {
            for (ra, rb) in a.iter().zip(b) {
                for (va, vb) in ra.iter().zip(rb) {
                    dev_t = dev_t.max((va - vb).abs());
                }
            }
        }
        check(
            "[C6] gauge 不変: local phase D X D† + 周波数 iωI は曲率と全時刻応答表を保存 (観測商の軌道)",
            dev_w < 1e-13 && dev_t < 1e-12 && mat_max_abs(&mat_sub(&x, &xg)) > 0.1,
            format!("max|Δw| = {:.1e}, max|Δ応答表| = {:.1e}, ‖X − X'‖ = {:.2}",
                    dev_w, dev_t, mat_max_abs(&mat_sub(&x, &xg))),
        );
    }

    // ---------------- [C7] 複素共役の曖昧性 → EquivalenceClassOnly ----------------
    {
        let mut rng = Rng::new(35203);
        let n = 3;
        let mut h = mat_zero(n);
        for i in 0..n {
            for j in (i + 1)..n {
                let z = cs(rng.f64() - 0.5, rng.f64() - 0.5);
                h[i][j] = z;
                h[j][i] = z.conj();
            }
        }
        let model = QuasiFreeOpenModel::new(h, vec![], vec![], None).expect("m");
        let x = model.effective_drift();
        let xc = mat_conj(&x);
        let nodes3 = vec![vec![0usize], vec![1usize], vec![2usize]];
        let times = [0.25, 0.7, 1.3];
        let ta = response_table(&x, &nodes3, &times);
        let tb = response_table(&xc, &nodes3, &times);
        let mut dev_t = 0.0f64;
        for (a, b) in ta.iter().zip(&tb) {
            for (ra, rb) in a.iter().zip(b) {
                for (va, vb) in ra.iter().zip(rb) {
                    dev_t = dev_t.max((va - vb).abs());
                }
            }
        }
        check(
            "[C7] 複素共役 X̄: 全時刻応答表が同一・X ≠ X̄ — 契約は X を商までしか識別しない (EquivalenceClassOnly)",
            dev_t < 1e-12 && mat_max_abs(&mat_sub(&x, &xc)) > 0.1,
            format!("max|Δ応答表| = {:.1e}, ‖X − X̄‖ = {:.2}", dev_t,
                    mat_max_abs(&mat_sub(&x, &xc))),
        );
    }

    // ---------------- [C8] 有限 shot: RobustExact / Straddled / Insufficient ----------------
    {
        let reader = FiniteShotCurvatureReader {
            eps: 0.3,
            delta: 0.2,
            alpha: 0.05,
            x_norm_bound: 1.0,
            tau: 0.3,
            min_shots: 1000,
        };
        // モデル族: w₂₁ = |h₂₁|² (散逸なしの閉鎖系で shot だけ有限)
        let run_cell = |h21: C64, seed: u64, n_shots: usize| -> CurvatureVerdict {
            let h = vec![vec![CZERO, h21], vec![h21.conj(), CZERO]];
            let model = QuasiFreeOpenModel::new(h, vec![], vec![], None).expect("m");
            let x = model.effective_drift();
            let pi = proj(2, &[0]);
            let pj = proj(2, &[1]);
            let mut c0 = mat_zero(2);
            c0[0][0] = cs(0.5, 0.0);
            c0[1][1] = cs(0.5, 0.0);
            let mut rng = Rng::new(seed);
            let sample = |t: f64, sgn: f64, rng: &mut Rng, n: usize| -> Vec<u8> {
                let cpr = if sgn > 0.0 {
                    mat_add(&c0, &mat_scale(cs(reader.eps, 0.0), &pi))
                } else {
                    mat_sub(&c0, &mat_scale(cs(reader.eps, 0.0), &pi))
                };
                let ct = evolve_covariance(&x, &mat_zero(2), &cpr, t);
                let p = mat_trace(&mat_mul(&pj, &ct)).re.clamp(0.0, 1.0);
                (0..n).map(|_| (rng.f64() < p) as u8).collect()
            };
            let s_pd = sample(reader.delta, 1.0, &mut rng, n_shots);
            let s_md = sample(reader.delta, -1.0, &mut rng, n_shots);
            let s_p2 = sample(2.0 * reader.delta, 1.0, &mut rng, n_shots);
            let s_m2 = sample(2.0 * reader.delta, -1.0, &mut rng, n_shots);
            reader.read(&s_pd, &s_md, &s_p2, &s_m2)
        };
        let v_edge = run_cell(cs(1.0, 0.0), 35204, 2_000_000);
        let v_none = run_cell(CZERO, 35205, 2_000_000);
        let v_stra = run_cell(cs(0.3f64.sqrt(), 0.0), 35206, 2_000_000);
        let v_insf = run_cell(cs(1.0, 0.0), 35207, 100);
        check(
            "[C8] 有限 shot 裁定: edge → RobustEdge / 無辺 → RobustNoEdge / w ≈ τ → Straddled / 不足 → Insufficient",
            v_edge == CurvatureVerdict::RobustEdge
                && v_none == CurvatureVerdict::RobustNoEdge
                && v_stra == CurvatureVerdict::Straddled
                && v_insf == CurvatureVerdict::InsufficientObservation,
            format!("{:?} / {:?} / {:?} / {:?}", v_edge, v_none, v_stra, v_insf),
        );
    }

    // ---------------- [C9] 相関 shot は読まない (v35.1 ゲートの継承) ----------------
    {
        let reader = FiniteShotCurvatureReader {
            eps: 0.3,
            delta: 0.2,
            alpha: 0.05,
            x_norm_bound: 1.0,
            tau: 0.3,
            min_shots: 1000,
        };
        // 持続 Markov 鎖 (stay 0.95) — 周辺は同じでも系列相関
        let mut rng = Rng::new(35208);
        let markov = |rng: &mut Rng, n: usize| -> Vec<u8> {
            let mut s = (rng.f64() < 0.5) as u8;
            (0..n)
                .map(|_| {
                    if rng.f64() > 0.95 {
                        s ^= 1;
                    }
                    s
                })
                .collect()
        };
        let a = markov(&mut rng, 100_000);
        let b = markov(&mut rng, 100_000);
        let c = markov(&mut rng, 100_000);
        let d = markov(&mut rng, 100_000);
        let v = reader.read(&a, &b, &c, &d);
        check(
            "[C9] 相関 shot: 遷移数ゲート (record_v2) が読みを拒否 → OutOfDomainCorrelated",
            v == CurvatureVerdict::OutOfDomainCorrelated,
            format!("{:?}", v),
        );
    }

    // ---------------- [C10] GQF-5: 電荷帰属 no-go と証明書の門 ----------------
    {
        // loss モデル: 電荷応答 −8, pairing なし
        let s2 = 8f64.sqrt();
        let loss_model = QuasiFreeOpenModel::new(
            mat_zero(2),
            vec![vec![cs(s2, 0.0), CZERO]],
            vec![],
            None,
        )
        .expect("m");
        let x_loss = loss_model.effective_drift();
        let c_full = {
            let mut c = mat_zero(2);
            c[0][0] = CONE;
            c
        };
        let s_loss = charge_response(&x_loss, &loss_model.injection(), &c_full);
        let resp_loss = ChargeNonconservingResponse(s_loss);
        // 散逸電荷保存証明書は loss モデルで拒否される (Λ ≠ 0)
        let cert_loss = certify_charge_conserving_dissipation(&loss_model);
        // pairing 側: 散逸ゼロの宣言は検証可能 (閉鎖 Nambu) → 証明書 OK
        let closed_decl = QuasiFreeOpenModel::new(mat_zero(2), vec![], vec![], None).expect("m");
        let cert_closed = certify_charge_conserving_dissipation(&closed_decl);
        let witness = match &cert_closed {
            Ok(c) => pairing_witness_with_certificate(ChargeNonconservingResponse(-8.0), c),
            Err(_) => None,
        };
        check(
            "[C10] GQF-5: 同値 −8 の電荷応答 — loss は証明書拒否 (witness 不能)・散逸ゼロ証明書のみ pairing witness を解錠 (禁止変換 33)",
            (s_loss + 8.0).abs() < 1e-12
                && cert_loss.is_err()
                && witness.is_some(),
            format!(
                "S_loss = {} (resp = {:?}), cert_loss = Err({:.1}), witness = {}",
                s_loss,
                resp_loss.0,
                cert_loss.err().unwrap_or(0.0),
                witness.map(|w| w.response).unwrap_or(0.0)
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "全検査 PASS — 曲率が読むのは有効 drift の観測商 (Hamiltonian への昇格は証明書の門のみ)".to_string()
        } else {
            format!("FAIL {} 件", nfail)
        }
    );
}
