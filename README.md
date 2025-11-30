# Ecological Inference Supplementary

Rscript for the supplementary material for the paper “An Accurate, Fast, and Scalable Ecological Inference Algorithm for the R×C Case,” by P. Ubilla, D. Hermosilla, and C. Thraves. 
- Layout: `src/` holds scripts and shared helpers; `output/` stores intermediate results such as `ei_instances`, `simulated_instances`; `figures/` receives generated PDF plots.

## Scripts
- Generate EI datasets: `Rscript src/run_instances.R`.
- Simulated-instance plots: `Rscript src/plot_errors.R` (Probability MAE boxplots, default `figures/sim_prob_error_3.pdf`) and `Rscript src/plot_time.R` (runtime vs candidates, default `figures/sim_runtime_3.pdf`) using `output/simulated_instances`.
- EI dataset wins by district: `Rscript src/plot_districts.R --out=figures/nz_relative_wins.pdf [--workers=8] [--limit-inst=N] [--parallel=true]`.
- Pairwise method pies: `Rscript src/pie_chart.R [outfile] [--field=EI_V] [--inst-like=ei_] [--limit-inst=N] [--workers=M] [--parallel=true] [--higher-better=true]` (default outfile `figures/pairwise_pie_grid.pdf`).
- Differences and p-values: `Rscript src/p_values.R <target_method> [field_to_test] [--inst-like=ei_] [--limit-inst=N] [--workers=M] [--parallel=true]`.
- EI instance table: `Rscript src/table.R <field_to_average> [outfile.tex] [workers] [--inst-like=ei_] [--limit-inst=N] [--parallel=true]` (also prints ASCII and writes XLSX next to the TEX file).
- Simulated-instance table: `Rscript src/table_simulated.R <field_to_average> [outfile.tex] [workers] [--inst-like=I] [--limit-inst=N] [--parallel=true]`.
- Scalability benchmark: `Rscript src/scalability.R [--parallel]` (writes `output/scalability_results*.csv`).
- EM start sensitivity demo: `Rscript src/starting_point.R` (runs a small grid of random starts and prints Frobenius norms).
