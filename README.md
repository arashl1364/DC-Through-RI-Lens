
# Replication files for "Discrete Choice in Marketing through the Lens of Rational Inattention"

Structured as follows:

## 1: Hierarchical Simulation (chapter 3)
Reproduces estimation results in Table 2 and 3.

## 2: Cost Identification (chapter 3)
Reproduces estimation results for Table 4, Figure 1 and Appendix C.

## 3: RI DCM Features (chapter 4)
Ordered by (sub) section; Reproduces figures and tables from chapter 4.

## Notes:
-For 1/2: Run "_main.R" files. Files automatically source necessary functions. ".RData" files contain outputs reported in section 3, replication with seed(66). <br>
-For 3  : Each subsection per folder. Corresponding figures/tables are organized in scripts. Needed functions to execute scripts are contained outside of subsection folders but are sourced within each figure/table script.

Packages needed:<br>
-bayesm<br>
-ggplot2<br>
-entropy<br>
-MASS<br>
-Rcpp<br>
-RcppArmadillo<br>

