/-
v6.8 アノマリー探索の機械検証 — 制限付き計算定理を Lean 4 の定理にする

対象は v3.1 の主張 (claims.yml: QRN-GAUGE-003, C2):
  「探索領域 (SU(3) 表現 ≤ 3/3̄, SU(2) ≤ 2, |6Y| ≤ 9, 多重項 ≤ 5, 成分 ≤ 15) 内で、
    5 アノマリー + Witten + カイラル性 + 全因子帯電を満たす集合は、
    標準模型 1 世代の同値軌道 (共役 × U(1) 反転の 4 通り) にちょうど一致する」

証明の構造:
  - 定理 sm_is_solution:      SM が全条件を満たす (カーネルの decide — 信頼基盤は Lean 核のみ)
  - 定理 sm_breaks_if_shifted: 超電荷を 1 目盛ずらすと落ちる (同上)
  - 定理 sm_minimal_unique:    領域の全列挙で「解は 4 個 = SM 軌道のみ」(native_decide)
    * native_decide は Lean コンパイラを信頼基盤に含む (核のみの decide では列挙が遅すぎる)。
      この選択は docs/uft-v6.8.md に明記する。
    * 列挙は燃料つき再帰で、燃料切れは結果を毒する (ok=false) — 定理が成立するなら
      燃料は尽きておらず、列挙は完全である (黙った打ち切りは定理を偽にする)。

実行: ~/.elan/bin/lean proofs/Anomaly.lean   (終了コード 0 = 全定理が検証済み)
-/

-- ---------------- 表現テーブル (v6.2 の Rust 実装と同一の整数データ) ----------------
-- id: 0:(1,1) 1:(1,2) 2:(3,1) 3:(3̄,1) 4:(3,2) 5:(3̄,2)
structure Rep where
  cd : Int -- 色次元
  wd : Int -- 弱次元
  a3 : Int -- SU(3)³ 係数 (SU(2) 成分あたり): 3 → +1, 3̄ → −1
  t3 : Int -- 2T(色): 3/3̄ → 1
  t2 : Int -- 2T(弱): 2 → 1
  cj : Nat -- 共役タイプ id

def reps : Array Rep := #[
  ⟨1, 1, 0, 0, 0, 0⟩,
  ⟨1, 2, 0, 0, 1, 1⟩,
  ⟨3, 1, 1, 1, 0, 3⟩,
  ⟨3, 1, -1, 1, 0, 2⟩,
  ⟨3, 2, 1, 1, 1, 5⟩,
  ⟨3, 2, -1, 1, 1, 4⟩]

def rep (t : Nat) : Rep := reps.getD t ⟨1, 1, 0, 0, 0, 0⟩

/-- 多重項 = (表現タイプ id, 超電荷 6Y)。-/
abbrev Mult := Nat × Int

def conjM (m : Mult) : Mult := ((rep m.1).cj, -m.2)
def yflipM (m : Mult) : Mult := (m.1, -m.2)
def compsOfT (t : Nat) : Int := (rep t).cd * (rep t).wd

-- ---------------- 無矛盾条件 (v3.1/v6.2 と同一) ----------------

/-- 質量項を許す対がない: 自己共役な多重項がなく、互いに共役な対もない。-/
def chiralOK : List Mult → Bool
  | [] => true
  | m :: rest =>
    conjM m != m && rest.all (fun x => x != conjM m) && chiralOK rest

/-- 5 アノマリー + Witten + 全因子帯電 + カイラル性。-/
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
    if r.wd == 2 then wit := wit + r.cd
    if r.cd > 1 then hasC := true
    if r.wd > 1 then hasW := true
    if n != 0 then hasY := true
  return su3cub == 0 && su3sq == 0 && su2sq == 0 && grav == 0 && cubic == 0
    && wit % 2 == 0 && hasC && hasW && hasY && chiralOK s

-- ---------------- 標準模型 1 世代とその同値軌道 ----------------

/-- SM 1 世代: e^c=(1,1)₆, L=(1,2)₋₃, u^c=(3̄,1)₋₄, d^c=(3̄,1)₂, Q=(3,2)₁ -/
def sm : List Mult := [(0, 6), (1, -3), (3, -4), (3, 2), (4, 1)]

def leM (a b : Mult) : Bool := a.1 < b.1 || (a.1 == b.1 && a.2 ≤ b.2)
def insertM (m : Mult) : List Mult → List Mult
  | [] => [m]
  | x :: xs => if leM m x then m :: x :: xs else x :: insertM m xs
def sortM (s : List Mult) : List Mult := s.foldr insertM []

/-- 物理的に同値な 4 通り (恒等 / 荷電共役 / U(1) 反転 / 両方)。
    |6Y|≤9 の領域では非自明な整数倍のスケーリングは範囲外なので軌道はこの 4 つ。-/
def smOrbit : List (List Mult) :=
  [sortM sm, sortM (sm.map conjM), sortM (sm.map yflipM),
   sortM (sm.map (fun m => yflipM (conjM m)))]

def inOrbit (s : List Mult) : Bool := smOrbit.contains (sortM s)

-- ---------------- 領域の全列挙 (v3.1 と同じ深さ優先・非減少原子列) ----------------

def natoms : Nat := 6 * 19
/-- 原子 i = (タイプ i/19, 超電荷 (i%19)−9)。非減少列として列挙すれば各多重集合を一度ずつ生成。-/
def atom (i : Nat) : Mult := (i / 19, (Int.ofNat (i % 19)) - 9)

/-- 深さ優先の全列挙 (native_decide の実行時間のため増分和で書く)。
    戻り値は「条件を満たす集合の数」だが、異常は数を**毒する**:
      - SM 軌道に入らない解を見つけたら +1,000,000
      - 燃料切れ (起きないはずだが) は +1,000,000,000
    したがって定理「走査結果 = 4」の成立は、(a) 解がちょうど 4 個、(b) 全てが SM 軌道、
    (c) 列挙が完全 (燃料は尽きていない)、の 3 つを同時に含意する。
    カイラル条件は単調 (共役対・自己共役は要素を足しても消えない) なので枝刈りに使う。
    増分和の各式は checkAll と同一であり、checkAll 自体は下の核検証定理が SM 点で照合する。-/
def walk :
    Nat → Nat → List Mult → Nat → Int → Int → Int → Int → Int → Int → Int →
    Bool → Bool → Bool → Nat
  | 0, _, _, _, _, _, _, _, _, _, _, _, _, _ => 1000000000 -- 燃料切れ = 毒
  | fuel + 1, i, acc, len, comps, a3, s3, s2, gr, cu, wit, hc, hw, hy =>
    if i < natoms then
      let m := atom i
      let r := rep m.1
      let c := r.cd * r.wd
      let take :=
        if len < 5 && comps + c ≤ 15
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
              if inOrbit acc' then 1 else 1000000 -- 軌道外の解 = 毒
            else 0
          hit + walk fuel i acc' (len + 1) (comps + c) a3' s3' s2' gr' cu' wit' hc' hw' hy'
        else 0
      take + walk fuel (i + 1) acc len comps a3 s3 s2 gr cu wit hc hw hy
    else 0

/-- v3.1 領域の走査結果。-/
def scanResult : Nat := walk 200 0 [] 0 0 0 0 0 0 0 0 false false false

-- ---------------- 定理 ----------------

/-- [核のみで検証] SM 1 世代は全ての無矛盾条件を満たす。-/
theorem sm_is_solution : checkAll sm = true := by decide

/-- [核のみで検証] 超電荷を 1 目盛ずらす (e^c: 6→5) と条件が破れる。-/
theorem sm_breaks_if_shifted :
    checkAll [(0, 5), (1, -3), (3, -4), (3, 2), (4, 1)] = false := by decide

/-- [native_decide で検証] v3.1 領域の全列挙: 条件を満たす集合はちょうど 4 個で、
    すべて SM の同値軌道に属する。すなわち同値を除いて SM 1 世代が唯一の解である。
    (4 個 = 恒等/共役/反転/両方の 4 表現。軌道外の解・燃料切れは走査結果を
     4 より大きく毒する設計なので、本定理の成立は列挙の完全性も含意する。) -/
theorem sm_minimal_unique : scanResult = 4 := by native_decide

-- 記録 (results/v68_lean.txt に残る)。scanResult の値は定理 sm_minimal_unique が
-- native_decide (コンパイル実行) で確定する — インタープリタでの再評価 (#eval) は
-- 桁違いに遅いため行わない。
#eval IO.println s!"SM 軌道 (4 通り) = {smOrbit}"
#eval IO.println "定理: sm_is_solution / sm_breaks_if_shifted (核 decide), sm_minimal_unique: scanResult = 4 (native_decide)"
#eval IO.println "このファイルがエラーなく通ったこと自体が全定理の検証である"
