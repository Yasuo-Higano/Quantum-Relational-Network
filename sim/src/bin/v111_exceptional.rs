//! v11.1 例外群走査の完結 — E₆ 以外の例外群のカイラル芯は空である
//!
//! v8.2 は E₆ の 27 のカイラル芯が SM 1 世代に一致することを示した (U(1) の階段の
//! 群論的傍証)。残高 4 の「例外群の直接走査」の残りは G₂, F₄, E₇, E₈ である。
//!
//! 鍵は表現論の確立された事実 (C0): **G₂, F₄, E₇, E₈ の全ての既約表現は自己共役**
//! (実または擬実) であり、自己共役表現のフェルミオンは常に vectorlike に対にできる
//! ため、カイラル芯 (共役対と完全中性を剥がした残り) は空になる。複素表現を持つ
//! 例外群は E₆ だけである。
//!
//! 本バイナリはこの事実を「引用」で済ませず、**低次表現のウェイト系を整数厳密に構成
//! して自己共役性 W = −W (多重集合) を機械検証する**:
//!   G₂: 7, 14 / F₄: 26, 52 / E₇: 56, 133 / E₈: 248
//! (それぞれ物質が入り得る最低次の表現と随伴。座標は 2 倍/整数化して i64 で厳密に扱う)
//!
//! 構成の正しさ自体も検査する:
//!   - 次元 (ウェイト数) が正しい
//!   - ルート系はノルム² が正しく 2 種以下 (単純レース群は 1 種)
//!   - 単純ルートによる Weyl 鏡映で閉じる (ルート系の公理)
//!   - ウェイト系も Weyl 鏡映で閉じる (表現の必要条件)
//! 陰性対照: SU(3) の 3 (複素表現) は W ≠ −W と判定されること。
//!
//! 結論: E₆ 以外の例外群はカイラル物質を作れない。v7.1 (2 本目 = B−L)・v8.2
//! (3 本目なし・E₆ 芯 = SM) と合わせ、「例外群の頂点」も E₆ で閉じる。

use uft_sim::*;

type V8 = [i64; 8];

fn neg(v: &V8) -> V8 {
    let mut o = [0i64; 8];
    for i in 0..8 {
        o[i] = -v[i];
    }
    o
}

fn dot(a: &V8, b: &V8) -> i64 {
    (0..8).map(|i| a[i] * b[i]).sum()
}

/// 多重集合として W = −W か (ソート比較; 厳密)
fn self_conjugate(w: &[V8]) -> bool {
    let mut a: Vec<V8> = w.to_vec();
    let mut b: Vec<V8> = w.iter().map(neg).collect();
    a.sort();
    b.sort();
    a == b
}

/// Weyl 鏡映 s_α(v) = v − 2(v·α)/(α·α) α — 整数座標で厳密 (割り切れることを検査)
fn reflect(v: &V8, alpha: &V8) -> Option<V8> {
    let num = 2 * dot(v, alpha);
    let den = dot(alpha, alpha);
    if num % den != 0 {
        return None;
    }
    let c = num / den;
    let mut o = *v;
    for i in 0..8 {
        o[i] -= c * alpha[i];
    }
    Some(o)
}

/// 集合が鏡映で閉じるか (多重集合として不変)
fn closed_under(w: &[V8], alphas: &[V8]) -> bool {
    let mut sorted: Vec<V8> = w.to_vec();
    sorted.sort();
    for a in alphas {
        let mut im: Vec<V8> = Vec::with_capacity(w.len());
        for v in w {
            match reflect(v, a) {
                Some(x) => im.push(x),
                None => return false,
            }
        }
        im.sort();
        if im != sorted {
            return false;
        }
    }
    true
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// 符号列挙 (mask のビット = 反転位置)
fn signs(base: &[i64], f: &mut dyn FnMut(&[i64])) {
    let n = base.len();
    for m in 0..(1usize << n) {
        let mut v: Vec<i64> = base.to_vec();
        for (i, x) in v.iter_mut().enumerate() {
            if (m >> i) & 1 == 1 {
                *x = -*x;
            }
        }
        f(&v);
    }
}

fn pad(v: &[i64]) -> V8 {
    let mut o = [0i64; 8];
    o[..v.len()].copy_from_slice(v);
    o
}

fn main() {
    self_test();
    println!("=== v11.1 例外群走査の完結: E₆ 以外のカイラル芯は空 ===\n");
    let mut all_ok = true;
    let check = |name: &str, ok: bool, all: &mut bool| {
        *all &= ok;
        println!("    {}  {}", name, pass(ok));
    };

    // ================= G₂ (3 次元和ゼロ座標) =================
    println!("[G₂] 7 (最低次) と 14 (随伴) — 座標: Σ=0 の 3 次元格子");
    let mut g2_short: Vec<V8> = Vec::new(); // 短ルート ±(1,−1,0) 巡回
    let mut g2_long: Vec<V8> = Vec::new(); // 長ルート ±(2,−1,−1) 巡回
    for (i, j) in [(0usize, 1usize), (1, 2), (0, 2)] {
        let mut v = [0i64; 8];
        v[i] = 1;
        v[j] = -1;
        g2_short.push(v);
        g2_short.push(neg(&v));
    }
    for p in [[2i64, -1, -1], [-1, 2, -1], [-1, -1, 2]] {
        let v = pad(&p);
        g2_long.push(v);
        g2_long.push(neg(&v));
    }
    let g2_roots: Vec<V8> = g2_short.iter().chain(g2_long.iter()).cloned().collect();
    let mut g2_7: Vec<V8> = g2_short.clone();
    g2_7.push([0i64; 8]); // 7 = 6 短ルート + 零
    let mut g2_14 = g2_roots.clone();
    g2_14.push([0i64; 8]);
    g2_14.push([0i64; 8]); // 14 = 12 ルート + 零 × rank 2
    // 単純ルート: α₁ = (1,−1,0) (短), α₂ = (−1,2,−1) (長)
    let g2_simple = vec![pad(&[1, -1, 0]), pad(&[-1, 2, -1])];
    check("次元: |7|=7, |14|=14", g2_7.len() == 7 && g2_14.len() == 14, &mut all_ok);
    check(
        "ルートのノルム²: 短 2 / 長 6 (比 3 = G₂ の刻印)",
        g2_short.iter().all(|r| dot(r, r) == 2) && g2_long.iter().all(|r| dot(r, r) == 6),
        &mut all_ok,
    );
    check(
        "Weyl 閉性: ルート系・7・14 とも鏡映で閉じる",
        closed_under(&g2_roots, &g2_simple)
            && closed_under(&g2_7, &g2_simple)
            && closed_under(&g2_14, &g2_simple),
        &mut all_ok,
    );
    check(
        "自己共役: 7 = 7̄, 14 = 14̄",
        self_conjugate(&g2_7) && self_conjugate(&g2_14),
        &mut all_ok,
    );

    // ================= F₄ (4 次元, 座標 2 倍で整数化) =================
    println!("\n[F₄] 26 (最低次) と 52 (随伴) — 座標: 2 倍格子");
    let mut f4_long: Vec<V8> = Vec::new(); // ±2e_i ± 2e_j (24 本, ノルム² 8)
    for i in 0..4usize {
        for j in (i + 1)..4 {
            for (si, sj) in [(1i64, 1i64), (1, -1), (-1, 1), (-1, -1)] {
                let mut v = [0i64; 8];
                v[i] = 2 * si;
                v[j] = 2 * sj;
                f4_long.push(v);
            }
        }
    }
    let mut f4_short: Vec<V8> = Vec::new(); // ±2e_i (8 本) + (±1,±1,±1,±1) (16 本) — ノルム² 4
    for i in 0..4usize {
        let mut v = [0i64; 8];
        v[i] = 2;
        f4_short.push(v);
        f4_short.push(neg(&v));
    }
    signs(&[1, 1, 1, 1], &mut |s| f4_short.push(pad(s)));
    let f4_roots: Vec<V8> = f4_long.iter().chain(f4_short.iter()).cloned().collect();
    let mut f4_26: Vec<V8> = f4_short.clone();
    f4_26.push([0i64; 8]);
    f4_26.push([0i64; 8]); // 26 = 24 短ルート + 零 × 2
    let mut f4_52 = f4_roots.clone();
    for _ in 0..4 {
        f4_52.push([0i64; 8]);
    }
    // 単純ルート (2 倍座標): α₁=2e₂−2e₃, α₂=2e₃−2e₄, α₃=2e₄, α₄=(e₁−e₂−e₃−e₄)
    let f4_simple = vec![
        pad(&[0, 2, -2, 0]),
        pad(&[0, 0, 2, -2]),
        pad(&[0, 0, 0, 2]),
        pad(&[1, -1, -1, -1]),
    ];
    check("次元: |26|=26, |52|=52", f4_26.len() == 26 && f4_52.len() == 52, &mut all_ok);
    check(
        "ルート数 48 (長 24 ノルム² 8 / 短 24 ノルム² 4)",
        f4_roots.len() == 48
            && f4_long.iter().all(|r| dot(r, r) == 8)
            && f4_short.iter().all(|r| dot(r, r) == 4),
        &mut all_ok,
    );
    check(
        "Weyl 閉性: ルート系・26・52 とも鏡映で閉じる",
        closed_under(&f4_roots, &f4_simple)
            && closed_under(&f4_26, &f4_simple)
            && closed_under(&f4_52, &f4_simple),
        &mut all_ok,
    );
    check(
        "自己共役: 26 = 26̄, 52 = 52̄",
        self_conjugate(&f4_26) && self_conjugate(&f4_52),
        &mut all_ok,
    );

    // ================= E₈ (8 次元, 座標 2 倍) =================
    println!("\n[E₈] 248 (随伴 = 最低次) — 座標: 2 倍格子");
    let mut e8_roots: Vec<V8> = Vec::new();
    for i in 0..8usize {
        for j in (i + 1)..8 {
            for (si, sj) in [(1i64, 1i64), (1, -1), (-1, 1), (-1, -1)] {
                let mut v = [0i64; 8];
                v[i] = 2 * si;
                v[j] = 2 * sj;
                e8_roots.push(v);
            }
        }
    }
    signs(&[1, 1, 1, 1, 1, 1, 1, 1], &mut |s| {
        let minus = s.iter().filter(|&&x| x < 0).count();
        if minus % 2 == 0 {
            e8_roots.push(pad(s));
        }
    });
    let mut e8_248 = e8_roots.clone();
    for _ in 0..8 {
        e8_248.push([0i64; 8]);
    }
    // 単純ルート (2 倍座標, Bourbaki): α₁=(1,−1,…,−1,1)/…: ここでは検査に十分な生成系として
    // 隣接差 2e_i−2e_{i+1} (i=1..6) と 2e₆+2e₇, 偶スピノル (1,1,1,1,1,-1,-1,-1)… の代わりに
    // **全ルートでの閉性**を直接検査する (生成系の選び方に依らない強い検査)
    check("ルート数 240, ノルム² 全て 8", e8_roots.len() == 240 && e8_roots.iter().all(|r| dot(r, r) == 8), &mut all_ok);
    check("次元: |248|=248", e8_248.len() == 248, &mut all_ok);
    check(
        "Weyl 閉性: 全 240 ルートの鏡映で 248 が閉じる",
        closed_under(&e8_248, &e8_roots),
        &mut all_ok,
    );
    check("自己共役: 248 = 248̄", self_conjugate(&e8_248), &mut all_ok);

    // ================= E₇ (E₈ 内で実現, 座標 2 倍) =================
    println!("\n[E₇] 56 (最低次) と 133 (随伴) — E₈ 内: α=2e₇+2e₈ に直交する部分系");
    let alpha: V8 = pad(&[0, 0, 0, 0, 0, 0, 2, 2]);
    let e7_roots: Vec<V8> = e8_roots
        .iter()
        .filter(|r| dot(r, &alpha) == 0)
        .cloned()
        .collect();
    // 56 の重み: {β − α/2 : β ∈ E₈ ルート, β·α = 4} — 2 倍座標では α/2 = (0..,1,1),
    // さらに全体を 2 倍して整数化: λ = 2β − α
    let e7_56: Vec<V8> = e8_roots
        .iter()
        .filter(|b| dot(b, &alpha) == 4)
        .map(|b| {
            let mut v = [0i64; 8];
            for i in 0..8 {
                v[i] = 2 * b[i] - alpha[i];
            }
            v
        })
        .collect();
    let mut e7_133 = e7_roots.clone();
    for _ in 0..7 {
        e7_133.push([0i64; 8]);
    }
    check("ルート数 126 (α 直交), |56|=56, |133|=133", e7_roots.len() == 126 && e7_56.len() == 56 && e7_133.len() == 133, &mut all_ok);
    check(
        "Weyl 閉性: E₇ の全ルート鏡映で 56 と 133 が閉じる",
        closed_under(&e7_56, &e7_roots) && closed_under(&e7_133, &e7_roots),
        &mut all_ok,
    );
    check(
        "自己共役: 56 = 56̄ (擬実), 133 = 133̄",
        self_conjugate(&e7_56) && self_conjugate(&e7_133),
        &mut all_ok,
    );

    // ================= 陰性対照: SU(3) の 3 は複素 =================
    println!("\n[対照] SU(3) の 3 (複素表現) — 同じ判定器が非自己共役を検出するか");
    let su3_3: Vec<V8> = vec![pad(&[2, -1, -1]), pad(&[-1, 2, -1]), pad(&[-1, -1, 2])];
    check(
        "3 ≠ 3̄ (判定器は盲目でない)",
        !self_conjugate(&su3_3),
        &mut all_ok,
    );
    // E₆ の複素性の代表 (v8.2 の 27 と整合): E₇ の 56 を E₆ に制限すると 27+27̄+1+1 —
    // ここでは C0 事実の引用に留める (v8.2 が 27 の分解を厳密に検査済み)

    // ================= 結論 =================
    println!("\n[結論]");
    println!("  検査した全表現 (G₂:7,14 / F₄:26,52 / E₇:56,133 / E₈:248) は W = −W (自己共役)。");
    println!("  自己共役表現のフェルミオンは vectorlike に対にできる (C0) ので、カイラル芯は空。");
    println!("  複素表現を持つ例外群は E₆ のみ (C0: 例外群で複素表現を持つのは E₆ だけ) であり、");
    println!("  v8.2 (E₆ の 27 のカイラル芯 = SM 1 世代) と合わせて:");
    println!("  *** 例外群の走査は閉じた — カイラル物質の入口は E₆ の 27、その芯は SM 1 世代 ***");

    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-GAUGE-015".into())),
        (
            "self_conjugate".into(),
            Json::Obj(vec![
                ("G2_7".into(), Json::Bool(self_conjugate(&g2_7))),
                ("G2_14".into(), Json::Bool(self_conjugate(&g2_14))),
                ("F4_26".into(), Json::Bool(self_conjugate(&f4_26))),
                ("F4_52".into(), Json::Bool(self_conjugate(&f4_52))),
                ("E7_56".into(), Json::Bool(self_conjugate(&e7_56))),
                ("E7_133".into(), Json::Bool(self_conjugate(&e7_133))),
                ("E8_248".into(), Json::Bool(self_conjugate(&e8_248))),
            ]),
        ),
        ("control_su3_complex".into(), Json::Bool(!self_conjugate(&su3_3))),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v111_exceptional.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {}", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
