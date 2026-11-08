/-
Copyright (c) 2026 NNS contributors. All rights reserved.
Released under GPL-3.0 license as described in the file LICENSE of OVVO-Financial/NNS.
Authors: NNS contributors
-/
module

public import NNSProofs.PartialMoments

/-!
# Co-partial moments

Formal model of `co_lpm`, `co_upm`, `d_lpm` and `d_upm` in
`src/partial_moments.cpp` for two equal-length samples and natural-number
degrees. The per-observation factors are exactly `lowerDev` / `upperDev`, so
ties sit on the lower side as in the C++:

| NNS     | `x` side          | `y` side          |
|---------|-------------------|-------------------|
| `CoLPM` | lower (`x ≤ tx`)  | lower (`y ≤ ty`)  |
| `CoUPM` | upper (`x > tx`)  | upper (`y > ty`)  |
| `DLPM`  | upper (`x > tx`)  | lower (`y ≤ ty`)  |
| `DUPM`  | lower (`x ≤ tx`)  | upper (`y > ty`)  |
-/

namespace NNS

open Finset

variable {n : ℕ}

@[expose] public noncomputable section

/-- Co-lower partial moment. -/
def CoLPM (dx dy : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, lowerDev dx tx (x i) * lowerDev dy ty (y i)) / n

/-- Co-upper partial moment. -/
def CoUPM (dx dy : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, upperDev dx tx (x i) * upperDev dy ty (y i)) / n

/-- Divergent lower partial moment: `x` above its target, `y` at or below. -/
def DLPM (dLower dUpper : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, upperDev dUpper tx (x i) * lowerDev dLower ty (y i)) / n

/-- Divergent upper partial moment: `x` at or below its target, `y` above. -/
def DUPM (dLower dUpper : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, lowerDev dLower tx (x i) * upperDev dUpper ty (y i)) / n

/-- Population co-moment about the targets: `(1/n) Σ (x - tx)(y - ty)`. -/
def coMoment (tx ty : ℝ) (x y : Fin n → ℝ) : ℝ :=
  (∑ i, (x i - tx) * (y i - ty)) / n

/-- Population covariance. -/
def covariance (x y : Fin n → ℝ) : ℝ := coMoment (mean x) (mean y) x y

lemma upperDev_sub_lowerDev_one (t v : ℝ) : upperDev 1 t v - lowerDev 1 t v = v - t := by
  have h := upperDev_add_neg_one_pow_mul_lowerDev 1 t v
  rw [pow_one, pow_one, neg_one_mul, ← sub_eq_add_neg] at h
  exact h

lemma lowerDev_zero_add_upperDev_zero (t v : ℝ) : lowerDev 0 t v + upperDev 0 t v = 1 := by
  have h := upperDev_add_neg_one_pow_mul_lowerDev 0 t v
  rw [pow_zero, one_mul, pow_zero, add_comm] at h
  exact h

/-- **Covariance decomposition.** At degree 1 the four co-partial moments
recover the co-moment about the targets:
`CoUPM + CoLPM - DLPM - DUPM = (1/n) Σ (x - tx)(y - ty)`. -/
theorem coMoment_decomposition (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoUPM 1 1 tx ty x y + CoLPM 1 1 tx ty x y - DLPM 1 1 tx ty x y - DUPM 1 1 tx ty x y
      = coMoment tx ty x y := by
  unfold CoUPM CoLPM DLPM DUPM coMoment
  rw [← add_div, ← sub_div, ← sub_div, ← sum_add_distrib, ← sum_sub_distrib,
    ← sum_sub_distrib]
  congr 1
  refine sum_congr rfl fun i _ => ?_
  rw [← upperDev_sub_lowerDev_one tx (x i), ← upperDev_sub_lowerDev_one ty (y i)]
  ring

/-- With the targets at the means, the decomposition gives the population
covariance. -/
theorem covariance_decomposition (x y : Fin n → ℝ) :
    CoUPM 1 1 (mean x) (mean y) x y + CoLPM 1 1 (mean x) (mean y) x y
      - DLPM 1 1 (mean x) (mean y) x y - DUPM 1 1 (mean x) (mean y) x y
      = covariance x y :=
  coMoment_decomposition (mean x) (mean y) x y

/-- **Quadrant partition.** At degree 0 the four co-partial moments are the
shares of observations in the four quadrants around `(tx, ty)`, so they sum
to `1`. -/
theorem quadrants_sum_eq_one (hn : 0 < n) (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoLPM 0 0 tx ty x y + CoUPM 0 0 tx ty x y + DLPM 0 0 tx ty x y + DUPM 0 0 tx ty x y
      = 1 := by
  unfold CoLPM CoUPM DLPM DUPM
  rw [← add_div, ← add_div, ← add_div, ← sum_add_distrib, ← sum_add_distrib,
    ← sum_add_distrib]
  have hterm : ∀ i : Fin n,
      lowerDev 0 tx (x i) * lowerDev 0 ty (y i) + upperDev 0 tx (x i) * upperDev 0 ty (y i)
        + upperDev 0 tx (x i) * lowerDev 0 ty (y i)
        + lowerDev 0 tx (x i) * upperDev 0 ty (y i) = 1 := fun i => by
    have hx := lowerDev_zero_add_upperDev_zero tx (x i)
    have hy := lowerDev_zero_add_upperDev_zero ty (y i)
    calc _ = (lowerDev 0 tx (x i) + upperDev 0 tx (x i))
              * (lowerDev 0 ty (y i) + upperDev 0 ty (y i)) := by ring
      _ = 1 := by rw [hx, hy, one_mul]
  rw [sum_congr rfl fun i _ => hterm i, sum_const, card_univ, Fintype.card_fin,
    nsmul_eq_mul, mul_one, div_self (Nat.cast_ne_zero.mpr hn.ne')]

/-- Co-partial moments of any degree are nonnegative. -/
theorem CoLPM_nonneg (dx dy : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    0 ≤ CoLPM dx dy tx ty x y :=
  div_nonneg (sum_nonneg fun _ _ =>
    mul_nonneg (lowerDev_nonneg _ _ _) (lowerDev_nonneg _ _ _)) n.cast_nonneg

theorem CoUPM_nonneg (dx dy : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    0 ≤ CoUPM dx dy tx ty x y :=
  div_nonneg (sum_nonneg fun _ _ =>
    mul_nonneg (upperDev_nonneg _ _ _) (upperDev_nonneg _ _ _)) n.cast_nonneg

/-- Symmetry: swapping the variables turns `DLPM` into `DUPM`. -/
theorem DLPM_swap (dl du : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    DLPM dl du tx ty x y = DUPM dl du ty tx y x := by
  unfold DLPM DUPM
  congr 1
  exact sum_congr rfl fun i _ => mul_comm _ _

/-- Co-moments are symmetric under swapping the variables. -/
theorem CoLPM_swap (dx dy : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoLPM dx dy tx ty x y = CoLPM dy dx ty tx y x := by
  unfold CoLPM
  congr 1
  exact sum_congr rfl fun i _ => mul_comm _ _

theorem CoUPM_swap (dx dy : ℕ) (tx ty : ℝ) (x y : Fin n → ℝ) :
    CoUPM dx dy tx ty x y = CoUPM dy dx ty tx y x := by
  unfold CoUPM
  congr 1
  exact sum_congr rfl fun i _ => mul_comm _ _

/-- A sample's co-lower partial moment with itself is a univariate LPM whose
degree is the sum of the two degrees. -/
theorem CoLPM_self (d e : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    CoLPM d e t t x x = LPM (d + e) t x := by
  unfold CoLPM LPM lowerDev
  congr 1
  refine sum_congr rfl fun i _ => ?_
  split_ifs <;> simp [pow_add]

/-- Upper analogue of `CoLPM_self`. -/
theorem CoUPM_self (d e : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    CoUPM d e t t x x = UPM (d + e) t x := by
  unfold CoUPM UPM upperDev
  congr 1
  refine sum_congr rfl fun i _ => ?_
  split_ifs <;> simp [pow_add]

/-- A sample never diverges from itself: `DLPM(x, x) = 0`. -/
theorem DLPM_self (dl du : ℕ) (t : ℝ) (x : Fin n → ℝ) : DLPM dl du t t x x = 0 := by
  unfold DLPM upperDev lowerDev
  rw [div_eq_zero_iff]
  left
  refine sum_eq_zero fun i _ => ?_
  split_ifs with h1 h2 <;> first | exact absurd h2 (not_le.mpr h1) | simp

end

end NNS
