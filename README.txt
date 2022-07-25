# ====================================================
# PRICING CLIMATE LINKERS
# ====================================================
# 
# Pauline Chikhani and Jean-Paul Renne
# 
# This version: July 2022.
# ====================================================


These codes replicate the results of the paper entitled "Pricing Climate Linkers," available at https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3881262



# ====================================================
A- How to run the codes?
# ====================================================
Run "main.R" to launch the estimation and produce the paper's charts & tables.

Before sourcing "main.R", modify the working directory (line 6 of main.R).

"main.R" launches calculates the model's solutions with initial parameters, and then recalculates the model with optimized parameters. 

If you want to re-optimize the parameters, switch "indic_estim" from 0 to 1.
 --- Weights of targets can be changed in "load_ini_model.R"

If you want to re-produce charts, switch "indic_plots_paper" from 0 to 1.
 --- Change the number of cores you want to dedicate to parallel computations on line 14 "nb.cores".
 --- Description of the plots is available from line 35 to 45.

If you want to re-produce tables, switch "indic_tables_paper" from 0 to 1.
 --- Description of the plots is available from line 52 to 55.

Note: using a laptop with a Apple M1 processor, CPU: 8 cores, launching the solutions of the optimized model with checks takes about 15 seconds, 
      one model ("model_solve") resolution 1 second.
Note: using a laptop with a Apple M1 processor, CPU: 8 cores, optimization of parameters takes about 30 minutes.
Note: using a laptop with a Apple M1 processor, CPU: 8 cores, producing all set of charts takes about 10 minutes.

# ====================================================
B- Organisation of codes
# ====================================================

The folders are as follows:

- "estimation" contains codes allowing for the estimation of the model (scripts of initial model + optimized model).
- "procedures" contains various procedures, organised in different scripts according to their use.
- "data" contains the data used by the charts. 
- "outputs" contains output charts and tables, as well as the scripts producing these output files.
