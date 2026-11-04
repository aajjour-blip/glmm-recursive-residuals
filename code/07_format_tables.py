"""Format manuscript Tables 11, B1, B2, B3 and supplementary Table S1 from saved outputs.

Run from the root of the supplement:  python3 code/07_format_tables.py

Inputs (written by the R scripts):
  data/ablation_epilepsy/partA_pvalues_by_variant.csv     (code/06_ablation_epilepsy.R)
  data/ablation_epilepsy/partA_component_pvalues.csv      (code/06_ablation_epilepsy.R)
  data/ablation_epilepsy/partB_rejection_rates.csv        (code/06_ablation_epilepsy.R)
  data/completed_simulations/Q1_Bernoulli/bernoulli_Q1_cell_summary.csv
                                                          (code/01_simulation_bernoulli_power.R)
  data/simulation_catalogue/simulation_catalogue_primary_B1000.csv (catalogue outputs)
Outputs:
  tables/table_ablation_epil_pcombo.tex               Table 11
  tables/table_ablation_epil_components.tex           Table B1
  tables/table_ablation_epil_components_refined.tex   Table B2
  tables/table_ablation_epil_rates.tex                Table B3
  tables/table_supplement_calibration.tex             Table S1
The script uses only the Python standard library and does not change any numerical value;
it formats the values computed in R.
"""
import collections
import csv
import os

B = 1000
ABL = os.path.join("data", "ablation_epilepsy")


def read_csv(path):
    """Return the rows of a CSV file as a list of dictionaries."""
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def fmt_p(p, calibrated=True):
    """Format a p-value: 0.001 for the smallest attainable calibrated value, <0.001 below."""
    if calibrated and abs(p - 1 / (B + 1)) < 1e-9:
        return "0.001"
    if p < 0.001:
        return "$<$0.001"
    return f"{p:.3f}"


def write(path, lines):
    """Write a list of lines to a file, terminated by a newline."""
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


pv = {(r["model"], r["key"]): float(r["p"])
      for r in read_csv(os.path.join(ABL, "partA_pvalues_by_variant.csv"))}
comp = collections.defaultdict(dict)
for r in read_csv(os.path.join(ABL, "partA_component_pvalues.csv")):
    comp[(r["model"], r["variant"])][r["statistic"]] = float(r["p"])
rates = {(r["scenario"], r["key"]): r
         for r in read_csv(os.path.join(ABL, "partB_rejection_rates.csv"))}

VARIANTS = [
    ("proposed", "Proposed", "Sequential residuals, refit reference, full set"),
    ("marginal", "(i) Residual", "Marginal residuals ($b_i$ not updated)"),
    ("fullpost", "", "Full-cluster posterior residuals"),
    ("plugin", "(ii) Reference", "Plug-in reference (no refit)"),
    ("known", "", "Known-parameter reference (i.i.d.\\ $N(0,1)$)"),
    ("neither", "(i) and (ii)", "Marginal residuals, plug-in reference"),
    ("indexed", "(iii) Statistics", "Indexed cumulative statistics only"),
    ("direct", "", "Direct distributional statistics only"),
    ("ks", "", "Kolmogorov--Smirnov test, uncalibrated"),
    ("sw", "", "Shapiro--Wilk test, uncalibrated"),
]
RULE_AFTER = {"proposed", "fullpost", "known", "neither"}
MODELS = ["Poisson_RI", "NB2_RI", "NB2_refined"]

# ---- Table 11: omnibus p-values of all variants -------------------------------------
lines = [
    r"\begin{table}[htbp]", r"\centering", r"\small",
    r"\caption{Component ablation on the epilepsy data. Each row replaces one component of"
    r" the proposed diagnostic and holds the others at the proposed choice. Entries are"
    r" calibrated omnibus $p$-values based on $B=1000$ simulated reference datasets, for"
    r" which 0.001 is the smallest attainable value; the last two rows are classical"
    r" uncalibrated $p$-values computed on the sequential residuals.}",
    r"\label{tab:ablation-epil}",
    r"\begin{tabular*}{\textwidth}{@{}ll@{\extracolsep{\fill}}ccc@{}}", r"\toprule",
    r"Component & Variant & Poisson RI & NB2 RI & NB2 refined \\", r"\midrule",
]
for key, block, label in VARIANTS:
    calibrated = key not in ("ks", "sw")
    cells = " & ".join(fmt_p(pv[(m, key)], calibrated) for m in MODELS)
    lines.append(f"{block} & {label} & {cells}" + r" \\")
    if key in RULE_AFTER:
        lines.append(r"\midrule")
lines += [r"\bottomrule", r"\end{tabular*}", r"\end{table}"]
write(os.path.join("tables", "table_ablation_epil_pcombo.tex"), lines)

# ---- Tables B1 and B2: component p-values ---------------------------------------------
INDEX_LABEL = {"period_num": "period", "pred_mean": "fitted mean", "base": "baseline count",
               "age": "age", "trt": "treatment"}
FIXED_LABEL = {"Tmax_r": r"$T_{\max r}$", "Tmax_C": r"$T_{\max C}$",
               "Tdisp": r"$T_{\mathrm{disp}}$", "Tlag": r"$T_{\mathrm{lag}}$"}
STAT_ORDER = ([f"{k}_{i}" for i in ["period_num", "pred_mean", "base", "age", "trt"]
               for k in ["Tinf", "T2"]] + ["Tmax_r", "Tmax_C", "Tdisp", "Tlag"])
COLS = ["proposed", "marginal", "fullpost", "plugin", "known"]


def stat_label(name):
    """LaTeX label of a component statistic name such as 'Tinf_base'."""
    if name in FIXED_LABEL:
        return FIXED_LABEL[name]
    kind, index = name.split("_", 1)
    symbol = r"$T_{\infty}$" if kind == "Tinf" else r"$T_{2}$"
    return f"{symbol} ({INDEX_LABEL[index]})"


def component_table(model, label, path, caption):
    """Write the component-level p-value table of one epilepsy model."""
    out = [r"\begin{table}[htbp]", r"\centering", r"\small", r"\caption{" + caption + "}",
           r"\label{" + label + "}",
           r"\begin{tabular*}{\textwidth}{@{}l@{\extracolsep{\fill}}ccccc@{}}", r"\toprule",
           r" & \multicolumn{3}{c}{Residual construction (refit reference)} &"
           r" \multicolumn{2}{c}{Reference (sequential residuals)} \\",
           r"\cmidrule(lr){2-4}\cmidrule(lr){5-6}",
           r"Statistic & Sequential & Marginal & Full posterior & Plug-in & Known-parameter \\",
           r"\midrule"]
    for name in STAT_ORDER:
        cells = " & ".join(fmt_p(comp[(model, c)][name]) for c in COLS)
        out.append(stat_label(name) + " & " + cells + r" \\")
    combo = " & ".join(fmt_p(pv[(model, c)]) for c in COLS)
    out += [r"\midrule", r"$p_{\mathrm{combo}}$ & " + combo + r" \\",
            r"\bottomrule", r"\end{tabular*}", r"\end{table}"]
    write(path, out)


COMMON = (r"under the three residual constructions (refit reference) and under the plug-in"
          r" and known-parameter references (sequential residuals); $B=1000$.")
component_table("NB2_RI", "tab:ablation-epil-components",
                os.path.join("tables", "table_ablation_epil_components.tex"),
                r"Component-level calibrated $p$-values for the NB2 random-intercept epilepsy"
                r" model " + COMMON)
component_table("NB2_refined", "tab:ablation-epil-components-refined",
                os.path.join("tables", "table_ablation_epil_components_refined.tex"),
                r"Component-level calibrated $p$-values for the refined NB2 epilepsy model "
                + COMMON + r" The statistics indexed by log baseline count, which are part of"
                r" the component set, coincide with those indexed by baseline count because the"
                r" two indices induce the same ordering.")

# ---- Table B3: design-anchored rejection rates -------------------------------------------
SCENARIOS = [("null", "Null"), ("mean", "Mean misspecification"),
             ("serial1", r"Serial, $\rho=0.25$"), ("serial2", r"Serial, $\rho=0.50$")]
lines = [
    r"\begin{table}[htbp]", r"\centering", r"\small",
    r"\caption{Rejection rates (\%) at the 5\% level, with Wilson 95\% intervals, for datasets"
    r" generated from the fitted epilepsy models on the observed design; 200 datasets per"
    r" scenario and $B=199$ reference replicates per dataset. Null: refined NB2 model"
    r" generated and fitted. Mean misspecification: refined NB2 model generated, NB2"
    r" random-intercept model fitted. Serial: NB2 random-intercept model with a lagged"
    r" log-count term generated, NB2 random-intercept model fitted.}",
    r"\label{tab:ablation-epil-rates}", r"\resizebox{\textwidth}{!}{%",
    r"\begin{tabular}{llcccc}", r"\toprule",
    "Component & Variant & " + " & ".join(h for _, h in SCENARIOS) + r" \\", r"\midrule",
]
for key, block, label in VARIANTS:
    cells = []
    for scen, _ in SCENARIOS:
        r = rates[(scen, key)]
        cells.append(f"{100 * float(r['rate']):.1f} ({100 * float(r['ci_low']):.1f},"
                     f" {100 * float(r['ci_high']):.1f})")
    lines.append(f"{block} & {label} & " + " & ".join(cells) + r" \\")
    if key in RULE_AFTER:
        lines.append(r"\midrule")
lines += [r"\bottomrule", r"\end{tabular}}", r"\end{table}"]
write(os.path.join("tables", "table_ablation_epil_rates.tex"), lines)

# ---- Table S1: calibration evidence ---------------------------------------------------
cells = [r for r in read_csv(os.path.join("data", "completed_simulations", "Q1_Bernoulli",
                                          "bernoulli_Q1_cell_summary.csv"))
         if r["row_type"] == "null_level"]
cells.sort(key=lambda r: (int(r["m"]), int(r["n_i"])))
catalogue = {r["scenario"]: r for r in read_csv(
    os.path.join("data", "simulation_catalogue", "simulation_catalogue_primary_B1000.csv"))}
n_total = sum(int(r["N_success"]) for r in cells)
x_total = sum(int(r["detected_n"]) for r in cells)
pooled = collections.OrderedDict()
for r in read_csv(os.path.join("data", "completed_simulations", "Q1_Bernoulli",
                               "Q1_summary_by_gamma_clean.csv")):
    if float(r["gamma"]) == 0:
        pooled = r
assert int(pooled["total_success"]) == n_total and int(pooled["total_detected"]) == x_total
lines = [
    r"\begin{table}[htbp]", r"\centering", r"\small",
    r"\caption{Supplementary calibration evidence: empirical null rejection rates under"
    r" correctly specified generalized linear mixed models. For the Bernoulli family, 500 null"
    r" datasets were generated per design cell; $N_{\mathrm{null}}$ is the number of successful"
    r" runs, and intervals are Wilson 95\% intervals. For Poisson and Gamma, the reported value"
    r" is the calibrated omnibus $p_{\mathrm{combo}}$ from the single correctly specified"
    r" scenario of the simulation catalogue ($B=1000$ refit replicates each). The pooled"
    r" Bernoulli rate is slightly above the nominal $5\%$ level; the intervals of two cells"
    r" with $m=50$ lie above $5\%$.}",
    r"\label{tab:supp-calibration}", r"\resizebox{\textwidth}{!}{%",
    r"\begin{tabular}{llrrrr}", r"\toprule",
    r"Family & Design $(m,\, n_i)$ & $N_{\mathrm{null}}$ & Empirical rejection rate &"
    r" \multicolumn{2}{c}{Wilson 95\% CI} \\",
    r"\cmidrule(lr){5-6}", r" &  &  &  & Lower & Upper \\", r"\midrule",
]
for r in cells:
    lines.append(f"Bernoulli & $m={r['m']},\\; n_i={r['n_i']}$ & {r['N_success']} &"
                 f" {100 * float(r['detection_rate']):.1f}\\% &"
                 f" {100 * float(r['ci_low']):.1f}\\% & {100 * float(r['ci_high']):.1f}\\% \\\\")
lines += [r"\midrule",
          f"Bernoulli & Pooled ($9$ cells) & ${n_total}$ &"
          f" ${100 * float(pooled['pooled_detection_rate']):.1f}\\%$ &"
          f" ${100 * float(pooled['ci_low']):.1f}\\%$ & ${100 * float(pooled['ci_high']):.1f}\\%$"
          r" \\", r"\midrule"]
for fam, scen in [("Poisson", "poisson_correct"), ("Gamma  ", "gamma_correct")]:
    p = float(catalogue[scen]["p_combo_omnibus"])
    lines.append(f"{fam} & $m=100,\\; n_i=8$ & $1^{{\\dag}}$ &"
                 f" $p_{{\\mathrm{{combo}}}} = {p:.3f}$ & \\multicolumn{{2}}{{c}}{{---}} \\\\")
lines += [r"\bottomrule", r"\end{tabular}}",
          r"\TableNote{$^{\dag}$Single calibration run from the simulation catalogue; repeated"
          r" null runs were not computed for the Poisson and Gamma families. The pooled"
          r" Bernoulli null rejection rate is $6.2\%$ (Wilson 95\% CI $5.6\%$--$7.0\%$) over"
          r" " + str(n_total) + r" successful runs spanning the nine design cells.}",
          r"\end{table}"]
write(os.path.join("tables", "table_supplement_calibration.tex"), lines)
print("Wrote Tables 11, B1, B2, B3 and S1 to tables/")
