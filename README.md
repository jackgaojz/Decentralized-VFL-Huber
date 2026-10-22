CODE GUIDE

1. Simulation studies (simulation/)

Main-text figures:
- stepsize.R: Figure 2.
- gamma.R: Figure 3.
- p1p2tau-1.R, p1p2tau-2.R, p1p2tau-3.R: Figure 4.
- m.R, n.R, d_cv.R, d.R: Figure 5.
- compare.cv.R, compare.R: Figure 6.

Appendix sections:
- signal-to-noise.R: "Sensitivity analysis with respect to signal-to-noise ratio".
- topology_effect_ring.R: "Results under a ring graph".
- topology_effect_star.R: "Results under a star graph".

2. Real-data analyses (real data/)

- realdata_customer_*.R: "Customer personality analysis dataset".
- realdata_communities_*.R: "Communities and crime dataset".

In these filenames, "cv" denotes cross-validation, "rep" denotes repeated
experiments, and "contam" denotes experiments with training-response contamination.
The "cv" scripts guide tuning parameter selection, while fixed parameter values are used in subsequent repeated experiments to reduce computation time.

3. Usage notes

Required R packages: MASS, caret, ggplot2, ggpubr, and latex2exp.
Before running a script, update its setwd() path to your local working directory.
For the real-data analyses, place the corresponding input files in that directory:
- Customer dataset: marketing_campaign.csv.
  Source: https://www.kaggle.com/datasets/imakash3011/customer-personality-analysis
- Communities dataset: communities.data and communities.names.
  Source: https://archive.ics.uci.edu/dataset/183/communities+and+crime
