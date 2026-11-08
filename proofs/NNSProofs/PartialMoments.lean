/-
Copyright (c) 2026 NNS contributors. All rights reserved.
Released under GPL-3.0 license as described in the file LICENSE of OVVO-Financial/NNS.
Authors: NNS contributors
-/
module

public import Mathlib

/-!
# Univariate partial moments

Formal model of the univariate kernels `lpm` / `upm` in
`src/partial_moments.cpp` (ports of the R package's `LPM_C` / `UPM_C`),
restricted to natural-number degrees.

A sample of size `n` is a function `x : Fin n → ℝ`. The conventions match the
C++ exactly:

* `LPM` counts ties: observation `i` contributes when `target - x i ≥ 0`.
* `UPM` is strict: observation `i` contributes when `x i - target > 0`.
* Degree `0` contributes `1` per qualifying observation (Lean's `a ^ 0 = 1`
  agrees with the C++ `degree == 0` branch).
* For `n = 0` Lean gives `0` (`a / 0 = 0`) while the C++ `lpm` / `upm` give
  `NaN` (`0.0 / 0`). Theorems whose right side is not `0` assume `0 < n`.
-/

namespace NNS

open Finset

@[expose] public noncomputable section

variable {n : ℕ}

/-- Contribution of one observation `v` to the lower partial moment. -/
def lowerDev (d : ℕ) (t v : ℝ) : ℝ := if v ≤ t then (t - v) ^ d else 0

/-- Contribution of one observation `v` to the upper partial moment. -/
def upperDev (d : ℕ) (t v : ℝ) : ℝ := if t < v then (v - t) ^ d else 0

/-- Sample mean. -/
def mean (x : Fin n → ℝ) : ℝ := (∑ i, x i) / n

/-- Lower partial moment of degree `d` about `t`: `LPM(d, t, x)` in NNS. -/
def LPM (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  (∑ i, lowerDev d t (x i)) / n

/-- Upper partial moment of degree `d` about `t`: `UPM(d, t, x)` in NNS. -/
def UPM (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  (∑ i, upperDev d t (x i)) / n

/-- Population central moment of degree `d` about `t`. -/
def centralMoment (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : ℝ :=
  (∑ i, (x i - t) ^ d) / n

/-! ## Pointwise lemmas -/

lemma lowerDev_nonneg (d : ℕ) (t v : ℝ) : 0 ≤ lowerDev d t v := by
  unfold lowerDev
  split_ifs with h
  · exact pow_nonneg (sub_nonneg.mpr h) d
  · exact le_refl 0

lemma upperDev_nonneg (d : ℕ) (t v : ℝ) : 0 ≤ upperDev d t v := by
  unfold upperDev
  split_ifs with h
  · exact pow_nonneg (sub_nonneg.mpr h.le) d
  · exact le_refl 0

/-- Every observation is on exactly one side of the target, so the signed sum of
its two contributions recovers the full deviation `(v - t) ^ d`. -/
lemma upperDev_add_neg_one_pow_mul_lowerDev (d : ℕ) (t v : ℝ) :
    upperDev d t v + (-1) ^ d * lowerDev d t v = (v - t) ^ d := by
  unfold upperDev lowerDev
  by_cases h : t < v
  · have h' : ¬ v ≤ t := not_le.mpr h
    simp [h, h']
  · have h' : v ≤ t := not_lt.mp h
    simp only [h, h', ↓reduceIte, zero_add]
    rw [← mul_pow]
    ring

/-! ## Sample identities -/

/-- **Master identity.** `UPM(d) + (-1)^d · LPM(d)` is the `d`-th central moment
about the target. Every identity below is a special case. -/
theorem UPM_add_neg_one_pow_mul_LPM (d : ℕ) (t : ℝ) (x : Fin n → ℝ) :
    UPM d t x + (-1) ^ d * LPM d t x = centralMoment d t x := by
  unfold UPM LPM centralMoment
  rw [mul_div_assoc', ← add_div, mul_sum, ← sum_add_distrib]
  congr 1
  exact sum_congr rfl fun i _ => upperDev_add_neg_one_pow_mul_lowerDev d t (x i)

/-- Degree 0: `LPM(0, t, x) + UPM(0, t, x) = 1`. -/
theorem LPM_zero_add_UPM_zero (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    LPM 0 t x + UPM 0 t x = 1 := by
  have h := UPM_add_neg_one_pow_mul_LPM 0 t x
  simp only [pow_zero, one_mul, centralMoment, sum_const, card_univ,
    Fintype.card_fin, nsmul_eq_mul, mul_one] at h
  rw [add_comm, h, div_self (Nat.cast_ne_zero.mpr hn.ne')]

/-- Degree 1: `UPM(1, t, x) - LPM(1, t, x) = mean(x) - t`. -/
theorem UPM_one_sub_LPM_one (hn : 0 < n) (t : ℝ) (x : Fin n → ℝ) :
    UPM 1 t x - LPM 1 t x = mean x - t := by
  have h := UPM_add_neg_one_pow_mul_LPM 1 t x
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  simp only [pow_one, neg_one_mul, ← sub_eq_add_neg, centralMoment, sum_sub_distrib,
    sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul] at h
  rw [h, mean, sub_div, mul_div_cancel_left₀ t hn']

/-- Even degree: `LPM(d) + UPM(d)` is the `d`-th central moment. -/
theorem LPM_add_UPM_of_even {d : ℕ} (hd : Even d) (t : ℝ) (x : Fin n → ℝ) :
    LPM d t x + UPM d t x = centralMoment d t x := by
  have h := UPM_add_neg_one_pow_mul_LPM d t x
  rw [hd.neg_one_pow, one_mul] at h
  rw [add_comm, h]

/-- Degree 2 about the mean: `LPM(2, μ, x) + UPM(2, μ, x)` is the population
variance. -/
theorem LPM_two_add_UPM_two_mean (x : Fin n → ℝ) :
    LPM 2 (mean x) x + UPM 2 (mean x) x = centralMoment 2 (mean x) x :=
  LPM_add_UPM_of_even even_two (mean x) x

/-- Partial moments are nonnegative. -/
theorem LPM_nonneg (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : 0 ≤ LPM d t x :=
  div_nonneg (sum_nonneg fun i _ => lowerDev_nonneg d t (x i)) n.cast_nonneg

theorem UPM_nonneg (d : ℕ) (t : ℝ) (x : Fin n → ℝ) : 0 ≤ UPM d t x :=
  div_nonneg (sum_nonneg fun i _ => upperDev_nonneg d t (x i)) n.cast_nonneg

/-- Degree-0 LPM is the empirical CDF: the share of observations `≤ t`. -/
theorem LPM_zero_eq_ecdf (t : ℝ) (x : Fin n → ℝ) :
    LPM 0 t x = (#{i | x i ≤ t} : ℝ) / n := by
  unfold LPM lowerDev
  simp only [pow_zero]
  rw [sum_boole]

/-- Degree-0 LPM (the empirical CDF) is monotone in the target. -/
theorem LPM_zero_monotone (x : Fin n → ℝ) : Monotone fun t => LPM 0 t x := by
  intro s t hst
  simp only [LPM_zero_eq_ecdf]
  gcongr

/-- Degree-1 LPM is monotone in the target. -/
theorem LPM_one_monotone (x : Fin n → ℝ) : Monotone fun t => LPM 1 t x := by
  intro s t hst
  dsimp only
  unfold LPM lowerDev
  gcongr with i
  simp only [pow_one]
  split_ifs with h1 h2 h2
  · linarith
  · exact absurd (le_trans h1 hst) h2
  · linarith
  · exact le_refl 0

end

end NNS
