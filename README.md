# Codes for the paper "Data-Driven Sample Average Approximation with Covariate Information"

The included <a href = "https://github.com/rohitkannan/DD-SAA/blob/master/Data-Driven%20SAA%20with%20Covariate%20Information%20(R1).pdf" target="_blank">PDF</a> is the latest version of the paper.

Results in the paper were generated using Julia 0.6.4, Python 3.8.5, JuMP 0.18.5, Gurobi 8.1.0, GLMNet 0.3.0, and NumPy 1.19.2. Most of the code should readily port over to later versions of Julia.


**Folder Contents:**
* The main files for the different methods lie in the "main_files" folder
* Implementations of the different solution methods are in the "solution_methods" folder
* Parameter settings for the different cases are in the "parameter_settings" folder
* The different instances for the cases are in the "instances" folder
* Shell scripts and submit files for HTCondor are in the "shell_scripts" folder
* Python scripts for plotting the results are in the "plotting_scripts" folder


**Case Details:**
* **SAA** ("FI-SAA")
  * case1, case4, and case5: sigma = 5 and omega = 1, omega = 2, and omega = 3
  * case2 and case3: omega = 1 and sigma = 10, sigma = 20
* **ER-SAA**
  * case1, case2, and case3: sigma = 5, omega = 1, and OLS regression without heteroscedasticity estimation for d_x = 3, d_x = 10, and d_x = 100
  * case4, case5, and case6: sigma = 5, omega = 1, and kNN regression without heteroscedasticity estimation for d_x = 3, d_x = 10, and d_x = 100
  * case7 and case8: sigma = 5, omega = 1, and Lasso regression without heteroscedasticity estimation for d_x = 10 and d_x = 100
  * case9 and case10: d_x = 10, omega = 1, and OLS regression without heteroscedasticity estimation for sigma = 10 and sigma = 20
  * case11 and case12: d_x = 10, omega = 1, and kNN regression without heteroscedasticity estimation for sigma = 10 and sigma = 20
  * case13, case15, and case17: d_x = 10, sigma = 5, and OLS regression with heteroscedasticity estimation for omega = 1, omega = 2, and omega = 3
  * case14 and case16: d_x = 10, sigma = 5, and OLS regression without heteroscedasticity estimation for omega = 2 and omega = 3
  * case18, case20, and case22: d_x = 10, sigma = 5, and kNN regression with heteroscedasticity estimation for omega = 1, omega = 2, and omega = 3
  * case19 and case21: d_x = 10, sigma = 5, and kNN regression without heteroscedasticity estimation for omega = 2 and omega = 3
* **kNN-SAA** (Bertsimas and Kallus)
  * case1, case2, and case3: sigma = 5 and omega = 1 for d_x = 3, d_x = 10, and d_x = 100
  * case4 and case5: d_x = 10 and omega = 1 for sigma = 10 and sigma = 20
  * case6 and case7: d_x = 10 and sigma = 5 for omega = 2 and omega = 3
* **J-SAA and J+-SAA**
  * case1 and case2: sigma = 5, omega = 1, and OLS regression without heteroscedasticity estimation for d_x = 10 and d_x = 100
* **Point Prediction**
  * case1, case2, and case3: sigma = 5, omega = 1, and OLS regression for d_x = 3, d_x = 10, and d_x = 100
  * case4 and case5: sigma = 5, omega = 1, and Lasso regression for d_x = 10 and d_x = 100
  * case6 and case7: d_x = 10, omega = 1, and OLS regression for sigma = 10 and sigma = 20
  * case8 and case9: d_x = 10, sigma = 5, and OLS regression for omega = 2 and omega = 3
* **N-SAA** ("naive SAA")
  * case1, case2, and case3: sigma = 5 and omega = 1, omega = 2, and omega = 3
