# mars_redox_2021

A collection of source code and scripts to reproduce the results of Wordsworth et al., 2021, Nature Geoscience.

PCM_LBL/ contains PCM_LBL, the line-by-line radiative-convective model used to calculation surface temperature as a function of total pressure, H2 content and solar luminosity. More details on the operation of this code are given in the README file in that directory.

The evol_model/ directory contains the matlab scripts to run the stochastic atmospheric evolution model that was used to produce Figs. 1-2 in the main paper and S2-S4 in the supplementary material. run_main.m is the main script.

steakley_impact_results/ contains the post-impact GCM climate modeling results created by Kathryn Steakley.

The .ac2 files are input scripts for The Geochemist's Workbench required to reproduce Figs. 3A and 4A in the main paper.

plot_TSR_figure.m plots Fig. 3B in the main paper.

plot_manganese_kinetics_figure.m plots Fig. 4B in the main paper.
