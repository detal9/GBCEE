# GBCEE
The files included in this folder should allow replication of all the simulations presented in the paper 
Talbot D, Beaudoin C (2022) A generalized double robust Bayesian model averaging approach to causal effect estimation with
application to the Study of Osteoporotic Fractures. Journal of Causal Inference.

The following files are packages or functions that were used to analyze the simulated data:
- BayesPen_1.0.tar.gz
- BCEE_1.3.1.tar.gz
- glider.R
- OAL_V3.R

The remaining files are used to generate and analyze the simulated data.
- scenario1.R, ..., scenario7.R are the files used to produce the results in Tables 2 and 6 of the article.
- scenario1_binary.R, scenario2_binary.R, scenario4_binary.R, scenario5_binary.R produce the results in Table 7
- computationalTime20191025.R produces the results mentioned in the last paragraph of Section 4.3
- positivity.R produces the results in Figure 1
- scenario1_supplemental.R, ..., scenario7_supplemental.R produce the results in Tables 13 and 14
- scenario2_binary_c.R, scenario2_cc.R, scenario2_exp.R produce the results in Tables 10, 8 & 9, and 11 & 12, respectively
- simulations_oralce.R produce the results in Appendix B
- scenario5_bootstrap produce the results mentioned in sixth paragraph of Section 4.3 concerning the use of bootstrap 
  in GBCEE.
  
Hopefully, this covers everything in the folder and the code is sufficiently well-documented to allow replication. Feel free
to contact me if you have any question.
