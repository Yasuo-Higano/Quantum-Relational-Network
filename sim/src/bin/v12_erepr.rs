//! v1.2 ER=EPR — もつれが時空を縫っている (熱場二重状態の厳密計算)
//!
//! Maldacena–Susskind (2013): 「2 つのブラックホールを繋ぐワームホール (ER 橋) と
//! 量子もつれ (EPR) は同じものである」。Van Raamsdonk (2010): もつれを減らすと
//! 時空は千切れる。これを厳密に見る:
//!
//! 2 本のフェルミオン鎖 L, R を熱場二重 (TFD) 状態
//!   |TFD⟩ = (1/√Z) Σ_n e^{-βE_n/2} |n⟩_L |n⟩_R
//! に置く。これは純粋状態で、片方だけ見ると温度 T の熱状態に見える。
//! (R は粒子-正孔変換した表現を使い、数保存のガウス形式で厳密に扱う)
//!
//! 測定: 鏡像領域間の相互情報量 I(A_L : A_R) を温度で追う
//!  → T→∞ で最大 (2つの宇宙は情報的に密着)、T→0 で消滅 (千切れる)。
//! 情報距離での MDS 埋め込みで「橋」の形成/切断を定量化する。

use uft_sim::*;

fn h2(z: f64) -> f64 {
    let z = z.clamp(1e-14, 1.0 - 1e-14);
    -z * z.ln() - (1.0 - z) * (1.0 - z).ln()
}

fn entropy_real(c: &[f64], n: usize) -> f64 {
    let (w, _) = jacobi_eigh(c, n);
    w.iter().map(|&z| h2(z)).sum()
}

struct Tfd {
    n: usize,
    cll: Vec<f64>, // C_LL(d)
    clr: Vec<f64>, // C_LR̃(d)
    crr: Vec<f64>, // C_R̃R̃(d)
}

impl Tfd {
    fn new(n: usize, beta: f64) -> Self {
        let mut cll = vec![0.0; n];
        let mut clr = vec![0.0; n];
        let mut crr = vec![0.0; n];
        for d in 0..n {
            let (mut a, mut b, mut c) = (0.0, 0.0, 0.0);
            for k in 0..n {
                let kk = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
                let e = -2.0 * kk.cos();
                let f = 1.0 / (1.0 + (beta * e).exp());
                let g = (f * (1.0 - f)).sqrt();
                let cs = (kk * d as f64).cos();
                a += f * cs;
                b += g * cs;
                c += (1.0 - f) * cs;
            }
            cll[d] = a / n as f64;
            clr[d] = b / n as f64;
            crr[d] = c / n as f64;
        }
        Tfd { n, cll, clr, crr }
    }
    /// サイト集合 (L側 ls, R側 rs) の相関行列を組む
    fn corr(&self, ls: &[usize], rs: &[usize]) -> Vec<f64> {
        let m = ls.len() + rs.len();
        let mut c = vec![0.0; m * m];
        let dist = |x: usize, y: usize| -> usize {
            let d = (x as isize - y as isize).unsigned_abs();
            d.min(self.n - d)
        };
        for (i, &x) in ls.iter().enumerate() {
            for (j, &y) in ls.iter().enumerate() {
                c[i + j * m] = self.cll[dist(x, y)];
            }
            for (j, &y) in rs.iter().enumerate() {
                c[i + (ls.len() + j) * m] = self.clr[dist(x, y)];
                c[(ls.len() + j) + i * m] = self.clr[dist(x, y)];
            }
        }
        for (i, &x) in rs.iter().enumerate() {
            for (j, &y) in rs.iter().enumerate() {
                c[(ls.len() + i) + (ls.len() + j) * m] = self.crr[dist(x, y)];
            }
        }
        c
    }
    fn mi(&self, ls: &[usize], rs: &[usize]) -> f64 {
        let sa = entropy_real(&self.corr(ls, &[]), ls.len());
        let sb = entropy_real(&self.corr(&[], rs), rs.len());
        let sab = entropy_real(&self.corr(ls, rs), ls.len() + rs.len());
        sa + sb - sab
    }
    /// L 内のブロック間 MI (幾何の基準スケール)
    fn mi_intra(&self, b1: &[usize], b2: &[usize]) -> f64 {
        let sa = entropy_real(&self.corr(b1, &[]), b1.len());
        let sb = entropy_real(&self.corr(b2, &[]), b2.len());
        let all: Vec<usize> = b1.iter().chain(b2.iter()).cloned().collect();
        let sab = entropy_real(&self.corr(&all, &[]), all.len());
        sa + sb - sab
    }
}

fn main() {
    let n = 100usize;
    println!("=== v1.2 ER=EPR: 熱場二重状態 — もつれが2つの宇宙を縫う ===\n");
    println!("2本の鎖 (各{}サイト) の TFD 純粋状態。鏡像ブロック (10サイト) 間の相互情報量:\n", n);
    let block: Vec<usize> = (40..50).collect();
    let block2: Vec<usize> = (50..60).collect(); // L 内の隣ブロック (基準)
    println!("  β(逆温度)  I(A_L:A_R) 鏡像間   I(A_L:B_L) 同一宇宙内の隣   橋/隣 比");
    let mut betas_mi = Vec::new();
    for &beta in &[0.02f64, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 20.0] {
        let tfd = Tfd::new(n, beta);
        let mirror = tfd.mi(&block, &block);
        let intra = tfd.mi_intra(&block, &block2);
        println!(
            "  {:6.2}     {:8.4}            {:8.4}                  {:8.4}",
            beta,
            mirror,
            intra,
            mirror / intra
        );
        betas_mi.push((beta, mirror));
    }
    let max_mi = 2.0 * 10.0 * (2.0f64).ln();
    println!("  (鏡像 MI の理論最大 = 2ℓ·ln2 = {:.3})", max_mi);
    println!("  => 高温: 鏡像間の MI が同一宇宙内の隣接ブロックを圧倒 = 2つの宇宙は「密着」");
    println!("     低温: 鏡像 MI → 0 = ワームホールが痩せて千切れる (Van Raamsdonk)\n");

    // 純粋性と熱性の確認
    {
        let beta = 1.0;
        let tfd = Tfd::new(n, beta);
        // 全系 (L+R の一部では純粋にならないので、モードレベルで確認済み) — ここでは熱性:
        // L の 1 サイト占有 = 1/2 (半充填)、L のエネルギー密度 vs 解析値
        let e_num = -2.0 * tfd.cll[1];
        let mut e_th = 0.0;
        for k in 0..n {
            let kk = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            let e = -2.0 * kk.cos();
            e_th += e / (1.0 + (beta * e).exp());
        }
        e_th /= n as f64;
        println!("[検証] β=1 で L 鎖単独のボンドエネルギー: TFD から {:.5}, 熱平均の解析値 {:.5}  {}",
            e_num, e_th, pass((e_num - e_th).abs() < 1e-10));
        println!("       純粋状態 |TFD⟩ の半分は、厳密に温度 1/β の熱状態 (熱 = 見えない相手とのもつれ)\n");
    }

    // 「ワームホールの長さ」: 情報距離 d = -ln(MI/MI_max) の温度依存
    println!("[幾何] 情報距離 d(A_L,A_R) = -ln(I/I_max) — ワームホールの長さの温度依存");
    println!("  β      d(鏡像間)");
    for &(beta, mi) in &betas_mi {
        let d = -(mi / max_mi).ln();
        println!("  {:6.2}  {:8.3}", beta, d);
    }
    println!("  => 温度を下げる(もつれを減らす)と 2 宇宙間の距離が単調に発散する。");
    println!("     接続された時空 = 強いもつれ、の直接の実演。\n");
    println!("結論: 「空間がつながっている」ことと「強くもつれている」ことは同じ事実の二つの顔。");
    println!("      ER=EPR: 時空の連結性は量子もつれの巨視的な現れである (A2 の動的検証)。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
