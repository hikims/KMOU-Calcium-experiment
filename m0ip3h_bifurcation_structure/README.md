# m0ip3h_bifurcation_structure — oscillation region, fixed points, Hopf boundary, relaxation regime

## Question
Where in parameter space does the stereotyped reference REF_ST oscillate, how is the cycle born (Hopf or a
hard onset with a coexisting stable fixed point), and how far from the Hopf boundary does the data window sit?
This is the "dynamical-state classification" the audit asked for (item 6-4) and the structural basis of the
stereotypy proposition (the cycle is a relaxation cycle: T · Im λ / 2π ≫ 1).

## Method
`m0ip3h_bifurcation_grid.c` (OpenMP over grid cells): per cell an RK4 run (dt 0.02; Tend raised to 25 ε_H
for slow cells; half-width interpolated; clipping counted) gives osc, T, amplitude, W50, duty, H range; the fixed
points are found by scanning g(C) = dC/dt at the explicit steady P*(C), H*(C) on a log grid (bisected), the
Jacobian by central differences, the eigenvalues from the cubic (bisection + deflation).  Class: 0 steady,
1 cycle + unstable FP, 2 cycle + stable FP (coexistence), 3 no cycle + unstable FP, +10 several FPs.
Grids: (S̄ ∈ [0.01, 100]) × (ε_H ∈ [10, 10⁴]) 41×41 at REF_ST; (β ∈ [0.5, 40]) × (q_off ∈ [0, 1]) 41×21 at S̄ = 1;
lines in S̄ (REF_ST and the previous reference), ε_H, β.

## Result
- (S̄, ε_H): one fixed point everywhere (no bistability of steady states); the cycle exists for ε_H above a
  boundary that rises from ≈ 10 at S̄ = 0.1 to ≈ 130 at S̄ ≥ 30.  For S̄ ≳ 0.3 the Hopf curve (Re λ = 0) lies
  above the oscillation boundary: a band of 70 cells has a cycle coexisting with a stable fixed point (hard onset).
  REF_ST (S̄ = 1, ε_H = 1000) is an order of magnitude above the boundary; T · Im λ/2π ≈ 7.5 (relaxation regime).
  At small S̄ the unstable eigenvalue is real (the rest state loses stability through a real eigenvalue and the
  cycle is not a Hopf cycle at all).  Period ≈ 0.08–0.5 ε_H; the amplitude also grows with ε_H (d).
- (β, q_off) at S̄ = 1: the cycle exists for 1.3 ≲ β ≲ 8 at every q_off (Hopf at β = 1.75 and 8.6 on the line);
  the duty is ≤ 0.1 over most of the region and the amplitude is largest at q_off = 0 (stereotyped setting).
- Along S̄: REF_ST period 517 → 57 (9×) over [10⁻³, 10³] with amplitude within ±12 %; the previous reference
  (q_off = 1) changes the amplitude 5×.  In the Azzoni dose window S̄ ∈ [0.10, 1.0] the period ratio is 2.3.
Numerics: no clipped states on either grid.

## Files
- `m0ip3h_bifurcation_grid.c` — `grid=1`, `line=1`, `point=1` (lists all fixed points and eigenvalues)
- `run_m0ip3h_bifurcation.sh` → `results/bif_grid_S_epsH_ref_st.csv`, `bif_grid_beta_qoff_S1.csv`,
  `bif_line_{S_ref_st,S_ref_prev,epsH_ref_st,beta_ref_st}.csv`
- `plot_m0ip3h_bifurcation_structure.py` → `results/bif_summary.csv`; `figures/m0ip3h_bifurcation_structure_main.*`
  (a: period map with oscillation and Hopf boundaries, b: duty map (β, q_off), c: lines in S̄),
  `figures/..._supp.*` (a: Re λ, b: relaxation index, c: amplitude (β, q_off), d: ε_H line, e: β line, f: classes)

## Run
    make && ./run_m0ip3h_bifurcation.sh && python3 plot_m0ip3h_bifurcation_structure.py     (~35 min on 2 threads)
