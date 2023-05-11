# bayessbw
The raw data is located in **data/all_days_df.csv**.

**SBC Code**

(1) **code/noss_data_simulation.R** contains the functions necessary to sample populations from the prior distribution.

(2) **code/noss_regniere_structured_reduced.cpp** contains the TMB code to run with tmbstan.

(3) **code/noss_sbc_run_reduced.R** runs the SBC analysis. Posterior distributions are written into **code/output**.


**Model Fitting and Weather Simulations**

(1) **code/noss_new_post.R** fits the model to the laboratory data.

(2) **code/noss_new_sims.R** pulls the weather data from the weathercan package and runs the weather simulations.
