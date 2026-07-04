/-
v7.5 アノマリー探索の機械検証 (2) — 大表現込みの v5.2 域を Lean 4 の定理にする

対象は v5.2/v6.2-R3 の主張 (claims.yml: QRN-GAUGE-007, C2):
  「拡張領域 (SU(3) 表現 ≤ 8, SU(2) ≤ 3, |6Y| ≤ 9, 多重項 ≤ 5, 成分 ≤ 16) 内で、
    無矛盾条件を満たすカイラル物質集合は、標準模型 1 世代の同値軌道 4 通りに限る」

構造は proofs/Anomaly.lean と同一 (毒値設計・燃料つき再帰・核検証の SM 定理 +
native_decide の全列挙定理)。表現テーブルが 6 → 10 種に拡大し、群論係数
A(6)=7, A(8)=0, 2T(6)=5, 2T(8)=6, 2T(3_w)=4 (Slansky) が加わる。
Witten 大域アノマリーに寄与するのは半整数アイソスピン (この域では二重項のみ —
三重項は整数アイソスピンなので寄与しない)。

実行: ~/.elan/bin/lean proofs/AnomalyBig.lean   (終了コード 0 = 全定理検証済み)
-/

structure Rep where
  cd : Int
  wd : Int
  a3 : Int -- SU(3)³ 係数 (SU(2) 成分あたり)
  t3 : Int -- 2T(色)
  t2 : Int -- 2T(弱)
  cj : Nat

-- id: 0:(1,1) 1:(1,2) 2:(1,3) 3:(3,1) 4:(3̄,1) 5:(3,2) 6:(3̄,2) 7:(6,1) 8:(6̄,1) 9:(8,1)
def reps : Array Rep := #[
  ⟨1, 1, 0, 0, 0, 0⟩,
  ⟨1, 2, 0, 0, 1, 1⟩,
  ⟨1, 3, 0, 0, 4, 2⟩,
  ⟨3, 1, 1, 1, 0, 4⟩,
  ⟨3, 1, -1, 1, 0, 3⟩,
  ⟨3, 2, 1, 1, 1, 6⟩,
  ⟨3, 2, -1, 1, 1, 5⟩,
  ⟨6, 1, 7, 5, 0, 8⟩,
  ⟨6, 1, -7, 5, 0, 7⟩,
  ⟨8, 1, 0, 6, 0, 9⟩]

def rep (t : Nat) : Rep := reps.getD t ⟨1, 1, 0, 0, 0, 0⟩

abbrev Mult := Nat × Int

def conjM (m : Mult) : Mult := ((rep m.1).cj, -m.2)
def yflipM (m : Mult) : Mult := (m.1, -m.2)
def compsOfT (t : Nat) : Int := (rep t).cd * (rep t).wd

def chiralOK : List Mult → Bool
  | [] => true
  | m :: rest =>
    conjM m != m && rest.all (fun x => x != conjM m) && chiralOK rest

/-- 5 アノマリー + Witten + 全因子帯電 + カイラル性 (v6.2 の Rust full_check と同一の式)。-/
def checkAll (s : List Mult) : Bool := Id.run do
  let mut su3cub : Int := 0
  let mut su3sq : Int := 0
  let mut su2sq : Int := 0
  let mut grav : Int := 0
  let mut cubic : Int := 0
  let mut wit : Int := 0
  let mut hasC := false
  let mut hasW := false
  let mut hasY := false
  for m in s do
    let r := rep m.1
    let n := m.2
    su3cub := su3cub + r.a3 * r.wd
    su3sq := su3sq + r.t3 * r.wd * n
    su2sq := su2sq + r.t2 * r.cd * n
    grav := grav + r.cd * r.wd * n
    cubic := cubic + r.cd * r.wd * n * n * n
    if r.wd == 2 then wit := wit + r.cd -- 半整数アイソスピン = 二重項のみ寄与
    if r.cd > 1 then hasC := true
    if r.wd > 1 then hasW := true
    if n != 0 then hasY := true
  return su3cub == 0 && su3sq == 0 && su2sq == 0 && grav == 0 && cubic == 0
    && wit % 2 == 0 && hasC && hasW && hasY && chiralOK s

/-- SM 1 世代 (10 種テーブルでの id): e^c=(1,1)₆, L=(1,2)₋₃, u^c=(3̄,1)₋₄, d^c=(3̄,1)₂, Q=(3,2)₁ -/
def sm : List Mult := [(0, 6), (1, -3), (4, -4), (4, 2), (5, 1)]

def leM (a b : Mult) : Bool := a.1 < b.1 || (a.1 == b.1 && a.2 ≤ b.2)
def insertM (m : Mult) : List Mult → List Mult
  | [] => [m]
  | x :: xs => if leM m x then m :: x :: xs else x :: insertM m xs
def sortM (s : List Mult) : List Mult := s.foldr insertM []

def smOrbit : List (List Mult) :=
  [sortM sm, sortM (sm.map conjM), sortM (sm.map yflipM),
   sortM (sm.map (fun m => yflipM (conjM m)))]

def inOrbit (s : List Mult) : Bool := smOrbit.contains (sortM s)

def natoms : Nat := 10 * 19
def atom (i : Nat) : Mult := (i / 19, (Int.ofNat (i % 19)) - 9)

/-- 深さ優先の全列挙 (増分和・カイラル単調枝刈り・毒値設計 — Anomaly.lean と同一)。
    軌道外の解 +10⁶、燃料切れ +10⁹ で毒するので、定理「= 4」の成立が
    (a) 解 4 個 (b) 全て SM 軌道 (c) 列挙の完全性 を同時に含意する。-/
def walk :
    Nat → Nat → List Mult → Nat → Int → Int → Int → Int → Int → Int → Int →
    Bool → Bool → Bool → Nat
  | 0, _, _, _, _, _, _, _, _, _, _, _, _, _ => 1000000000
  | fuel + 1, i, acc, len, comps, a3, s3, s2, gr, cu, wit, hc, hw, hy =>
    if i < natoms then
      let m := atom i
      let r := rep m.1
      let c := r.cd * r.wd
      let take :=
        if len < 5 && comps + c ≤ 16
            && conjM m != m && acc.all (fun x => x != conjM m) then
          let n := m.2
          let a3' := a3 + r.a3 * r.wd
          let s3' := s3 + r.t3 * r.wd * n
          let s2' := s2 + r.t2 * r.cd * n
          let gr' := gr + r.cd * r.wd * n
          let cu' := cu + r.cd * r.wd * n * n * n
          let wit' := wit + (if r.wd == 2 then r.cd else 0)
          let hc' := hc || r.cd > 1
          let hw' := hw || r.wd > 1
          let hy' := hy || n != 0
          let acc' := m :: acc
          let hit :=
            if a3' == 0 && s3' == 0 && s2' == 0 && gr' == 0 && cu' == 0
                && wit' % 2 == 0 && hc' && hw' && hy' then
              if inOrbit acc' then 1 else 1000000
            else 0
          hit + walk fuel i acc' (len + 1) (comps + c) a3' s3' s2' gr' cu' wit' hc' hw' hy'
        else 0
      take + walk fuel (i + 1) acc len comps a3 s3 s2 gr cu wit hc hw hy
    else 0

def scanResult : Nat := walk 250 0 [] 0 0 0 0 0 0 0 0 false false false

/-- [核のみで検証] SM は拡張域でも全条件を満たす。-/
theorem sm_is_solution : checkAll sm = true := by decide

/-- [核のみで検証] 大表現の反例: (6,1) を含む単純な組は SU(3)³ が消えない。-/
theorem sextet_pair_fails :
    checkAll [(7, 1), (4, -4), (4, 2), (1, -3), (0, 6)] = false := by decide

/-- [native_decide で検証] v5.2 域 (大表現込み・|6Y|≤9・≤5 多重項・≤16 成分) の全列挙:
    解はちょうど 4 個で、すべて SM の同値軌道に属する — 孤立性は表現拡大に頑健。-/
theorem sm_unique_with_big_reps : scanResult = 4 := by native_decide

#eval IO.println s!"SM 軌道 (4 通り) = {smOrbit}"
#eval IO.println "定理: sm_is_solution / sextet_pair_fails (核 decide), sm_unique_with_big_reps: scanResult = 4 (native_decide)"
#eval IO.println "このファイルがエラーなく通ったこと自体が全定理の検証である"
