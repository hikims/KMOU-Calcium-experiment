"""
fig4_revised.py — Fig. 5 of the manuscript (file fig4_main_first_passage_theory) without the points that the revision
withdraws: full-model points with nu > 1.2 and reduced Monte-Carlo points with nu > 1 were computed with the
Euler-Maruyama scheme with reflection, which is inaccurate when 0 is a regular boundary (e1, e1b).
Uses the panel functions of the original folders through assemble_paper_figures; nothing else is changed.
Output: figures/fig4_main_first_passage_theory.{png,pdf,eps}
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ASM = os.path.join(HERE, "..", "assemble_paper_figures"); sys.path.insert(0, ASM)
import assemble_paper_figures as A
A.FIG = os.path.join(HERE, "figures"); os.makedirs(A.FIG, exist_ok=True)
A._patch_multiprocessing()
path, ns = A.load_folder("m1_first_passage_theory")
print("full-model points:", len(ns["R"]), "->", int((ns["R"].nu <= 1.2).sum()), "; reduced Monte-Carlo points:", len(ns["mc"]), "->", int((ns["mc"].nu <= 1.0).sum()))
ns["R"] = ns["R"][ns["R"].nu <= 1.2].copy(); ns["mc"] = ns["mc"][ns["mc"].nu <= 1.0].copy()
spec = [s for s in A.MAIN_FIGS if s[0] == "fig4_main_first_passage_theory"][0]
print(A.render(spec))
