/-
v14.4 アノマリー探索の機械検証 (3): v4.3 域の配列ベース全列挙

対象は v4.3/v6.2 の主張 (claims.yml: QRN-GAUGE-006, C2):
  「探索領域 (小表現, |6Y| ≤ 9, 多重項 ≤ 6, 成分 ≤ 24) 内で、無矛盾カイラル集合の
    正準形の成分数分布はちょうど {15 ↦ 1, 16 ↦ 8, 24 ↦ 459} (計 468) である」
  — SM (15 成分) は最小・一意で、16 成分の隣人 8 つの先は 24 成分まで空白 (孤立)。

従来の燃料つき再帰 walk (Anomaly.lean 方式) ではこの領域は ~10¹¹ 呼び出し (~17 時間)
と見積もられていた。本ファイルは**配列ベースの中間一致 (MITM)** で書き直す:

  - 半分列挙: 非減少な原子列 (長さ ≤ 3) を for ループの直積で構築 —
    **完全性は Array.range の構造から自明** (燃料も毒も不要)
  - 結合: 5 アノマリー和を単一 Int にパックした鍵で HashMap 結合
    (A 側の鍵 + B 側の鍵 = 0、Witten パリティ一致、接合部の非減少、成分 ≤ 24)
  - 候補ごとに全条件 (カイラル性・全因子帯電) を再検査し、正準形 (共役 × U(1) 反転
    の 4 変換の最小) で重複除去してから成分数分布を数える

定理 v43_spectrum の成立は (a) 分布がちょうど {15:1, 16:8, 24:459}、
(b) それ以外の成分数の解が存在しない (poison 項が 0)、を同時に含意する。
列挙の完全性は構造的 (全ループが固定範囲) なので、燃料切れの毒は不要になった。

実行: ~/.elan/bin/lean proofs/AnomalyArray.lean   (終了コード 0 = 全定理が検証済み)
予測実行時間: 半分列挙 ~26 万 × 2 + 結合 — native_decide で数十秒〜数分。
-/

import Std.Data.HashMap

open Std

-- ---------------- 表現テーブル (Anomaly.lean と同一) ----------------
structure Rep where
  cd : Int
  wd : Int
  a3 : Int
  t3 : Int
  t2 : Int
  cj : Nat

def reps : Array Rep := #[
  ⟨1, 1, 0, 0, 0, 0⟩,
  ⟨1, 2, 0, 0, 1, 1⟩,
  ⟨3, 1, 1, 1, 0, 3⟩,
  ⟨3, 1, -1, 1, 0, 2⟩,
  ⟨3, 2, 1, 1, 1, 5⟩,
  ⟨3, 2, -1, 1, 1, 4⟩]

def rep (t : Nat) : Rep := reps.getD t ⟨1, 1, 0, 0, 0, 0⟩

/-- 原子 = (表現タイプ, 超電荷 6Y) をひとつの添字 0..113 で表す。 -/
def natoms : Nat := 6 * 19
def atomT (i : Nat) : Nat := i / 19
def atomY (i : Nat) : Int := (Int.ofNat (i % 19)) - 9

/-- 原子の共役: タイプは cj、超電荷は −y (添字で 18 − (i%19))。 -/
def conjA (i : Nat) : Nat := (rep (atomT i)).cj * 19 + (18 - i % 19)
/-- U(1) 反転: 超電荷のみ −y。 -/
def yflipA (i : Nat) : Nat := (atomT i) * 19 + (18 - i % 19)

def compsA (i : Nat) : Int := (rep (atomT i)).cd * (rep (atomT i)).wd

-- ---------------- 無矛盾条件 (原子添字のリスト上で; Anomaly.lean と同一の式) ----------------

def chiralOK : List Nat → Bool
  | [] => true
  | m :: rest =>
    conjA m != m && rest.all (fun x => x != conjA m) && chiralOK rest

def checkAll (s : List Nat) : Bool := Id.run do
  let mut su3cub : Int := 0
  let mut su3sq : Int := 0
  let mut su2sq : Int := 0
  let mut grav : Int := 0
  let mut cubic : Int := 0
  let mut wit : Int := 0
  let mut hasC := false
  let mut hasW := false
  let mut hasY := false
  for i in s do
    let r := rep (atomT i)
    let n := atomY i
    su3cub := su3cub + r.a3 * r.wd
    su3sq := su3sq + r.t3 * r.wd * n
    su2sq := su2sq + r.t2 * r.cd * n
    grav := grav + r.cd * r.wd * n
    cubic := cubic + r.cd * r.wd * n * n * n
    if r.wd == 2 then wit := wit + r.cd
    if r.cd > 1 then hasC := true
    if r.wd > 1 then hasW := true
    if n != 0 then hasY := true
  return su3cub == 0 && su3sq == 0 && su2sq == 0 && grav == 0 && cubic == 0
    && wit % 2 == 0 && hasC && hasW && hasY && chiralOK s

-- ---------------- 半分列挙 ----------------

/-- 半分 = (原子添字列 [非減少], 末尾添字, 成分和, パック鍵, Witten パリティ)。
    鍵は 5 つのアノマリー和を安全な基数でパックした単一 Int:
    範囲 (≤ 3 原子): |a3|≤6, |s3|≤54, |s2|≤54, |gr|≤162, |cu|≤13122。 -/
structure Half where
  idxs : List Nat
  last : Nat
  comps : Int
  key : Int
  witp : Int

def packKey (a3 s3 s2 gr cu : Int) : Int :=
  (((a3 * 200 + s3) * 200 + s2) * 500 + gr) * 30000 + cu

/-- 原子 1 つ分の (a3, s3, s2, gr, cu, wit) 寄与。 -/
def contrib (i : Nat) : Int × Int × Int × Int × Int × Int :=
  let r := rep (atomT i)
  let n := atomY i
  (r.a3 * r.wd, r.t3 * r.wd * n, r.t2 * r.cd * n, r.cd * r.wd * n,
   r.cd * r.wd * n * n * n, if r.wd == 2 then r.cd else 0)

/-- 長さちょうど len (1..3) の非減少原子列を全列挙 (成分和 ≤ 24)。
    for ループの直積なので完全性は構造的。 -/
def halves (len : Nat) : Array Half := Id.run do
  let mut out : Array Half := #[]
  if len == 1 then
    for i in [0:natoms] do
      let (a3, s3, s2, gr, cu, w) := contrib i
      out := out.push ⟨[i], i, compsA i, packKey a3 s3 s2 gr cu, w % 2⟩
  else if len == 2 then
    for i in [0:natoms] do
      for j in [i:natoms] do
        let c := compsA i + compsA j
        if c ≤ 24 then
          let (a1, b1, c1, d1, e1, w1) := contrib i
          let (a2, b2, c2, d2, e2, w2) := contrib j
          out := out.push ⟨[i, j], j, c,
            packKey (a1+a2) (b1+b2) (c1+c2) (d1+d2) (e1+e2), (w1+w2) % 2⟩
  else
    for i in [0:natoms] do
      for j in [i:natoms] do
        for k in [j:natoms] do
          let c := compsA i + compsA j + compsA k
          if c ≤ 24 then
            let (a1, b1, c1, d1, e1, w1) := contrib i
            let (a2, b2, c2, d2, e2, w2) := contrib j
            let (a3', b3, c3, d3, e3, w3) := contrib k
            out := out.push ⟨[i, j, k], k, c,
              packKey (a1+a2+a3') (b1+b2+b3) (c1+c2+c3) (d1+d2+d3) (e1+e2+e3),
              (w1+w2+w3) % 2⟩
  return out

-- ---------------- 正準形 (共役 × U(1) 反転の 4 変換の最小) ----------------

def leL : List Nat → List Nat → Bool
  | [], [] => true
  | [], _ => true
  | _, [] => false
  | a :: as, b :: bs => a < b || (a == b && leL as bs)

def insertN (m : Nat) : List Nat → List Nat
  | [] => [m]
  | x :: xs => if m ≤ x then m :: x :: xs else x :: insertN m xs
def sortN (s : List Nat) : List Nat := s.foldr insertN []

/-- U(1) 電荷の gcd 約分 (正準形の「整数倍で不変」— v6.2 の Rust 正準形と同一)。
    g = gcd(|y|) ≥ 2 なら全電荷を g で割る (割り切れは gcd の定義から厳密)。 -/
def gcdReduce (s : List Nat) : List Nat := Id.run do
  let mut g : Nat := 0
  for i in s do
    g := Nat.gcd g (atomY i).natAbs
  if g ≤ 1 then
    return s
  else
    return s.map (fun i => (atomT i) * 19 + ((atomY i) / (Int.ofNat g) + 9).toNat)

def canon (s : List Nat) : List Nat := Id.run do
  let r := gcdReduce s
  let cands := [sortN r, sortN (r.map conjA), sortN (r.map yflipA),
                sortN (r.map (fun i => yflipA (conjA i)))]
  let mut best := sortN r
  for c in cands do
    if leL c best && c != best then
      best := c
  return best

-- ---------------- MITM 走査 ----------------

/-- 走査結果: (15 成分の正準解数, 16 成分, 24 成分, その他の成分数の解数 [毒]) -/
def arrayScan : Nat × Nat × Nat × Nat := Id.run do
  -- 半分の事前計算
  let h1 := halves 1
  let h2 := halves 2
  let h3 := halves 3
  -- B 側 (鍵 → Half 配列) の HashMap
  let mkMap (hs : Array Half) : HashMap Int (Array Half) := Id.run do
    let mut m : HashMap Int (Array Half) := {}
    for h in hs do
      m := m.insert h.key ((m.getD h.key #[]).push h)
    return m
  let m1 := mkMap h1
  let m2 := mkMap h2
  let m3 := mkMap h3
  -- 解の収集 (正準形)
  let mut sols : Array (List Nat) := #[]
  -- 分割 (a,b): L=1 → (1,0), L=2 → (1,1), L=3 → (2,1), L=4 → (2,2), L=5 → (3,2), L=6 → (3,3)
  -- (a, 0): 単独で完全な解
  for hA in h1 do
    if checkAll hA.idxs then
      sols := sols.push (canon hA.idxs)
  -- (a, b) 結合の共通処理
  let join (has : Array Half) (mb : HashMap Int (Array Half)) :
      Array (List Nat) := Id.run do
    let mut found : Array (List Nat) := #[]
    for hA in has do
      match mb.get? (-hA.key) with
      | none => pure ()
      | some bs =>
        for hB in bs do
          -- 接合部の非減少 (B の先頭 ≥ A の末尾)・成分・パリティ
          let bFirst := hB.idxs.headD 0
          if hA.last ≤ bFirst && hA.comps + hB.comps ≤ 24
              && (hA.witp + hB.witp) % 2 == 0 then
            let full := hA.idxs ++ hB.idxs
            if checkAll full then
              found := found.push (canon full)
    return found
  sols := sols ++ join h1 m1 ++ join h2 m1 ++ join h2 m2 ++ join h3 m2 ++ join h3 m3
  -- 重複除去 (ソートして走査)
  let sorted := sols.qsort (fun a b => leL a b && a != b)
  let mut n15 := 0
  let mut n16 := 0
  let mut n24 := 0
  let mut other := 0
  let mut prev : List Nat := []
  let mut first := true
  for s in sorted do
    if first || s != prev then
      first := false
      prev := s
      let c : Int := s.foldl (fun acc i => acc + compsA i) 0
      if c == 15 then n15 := n15 + 1
      else if c == 16 then n16 := n16 + 1
      else if c == 24 then n24 := n24 + 1
      else other := other + 1
  return (n15, n16, n24, other)

-- ---------------- 定理 ----------------

/-- SM 1 世代 (原子添字): e^c=(0, y=6→添字15), L=(1,−3→添字25), u^c=(3̄=タイプ3,−4→添字62),
    d^c=(タイプ3, 2→添字68), Q=(タイプ4, 1→添字86)。 -/
def smIdx : List Nat := [15, 25, 62, 68, 86]

/-- [核のみで検証] 添字エンコードの検算: SM は全条件を満たす。 -/
theorem sm_ok_in_index_encoding : checkAll smIdx = true := by decide

/-- [native_decide で検証] v4.3 域の配列ベース全列挙:
    正準形解 (並べ替え・共役・U(1) 反転・整数倍で同一視) の成分数分布は
    ちょうど {15 ↦ 1, 16 ↦ 8, 24 ↦ 459}、他の成分数は 0。
    完全性は構造的 (全ループが固定範囲の直積) — 燃料の毒は不要になった。 -/
theorem v43_spectrum : arrayScan = (1, 8, 459, 0) := by native_decide

#eval IO.println "定理: sm_ok_in_index_encoding (核 decide), v43_spectrum: arrayScan = (1, 8, 459, 0) (native_decide)"
#eval IO.println "このファイルがエラーなく通ったこと自体が全定理の検証である"
