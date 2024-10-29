# UV wing pattern diversification in Papilionidae

This repository contains the R code to run statistical analysis on phenotypic coordinates of UV wing colour pattern of 77 Papilionidae species.

Within the code folder :

The file "brightness.R" runs the phylogenetically paired t-test to test for difference in average UV brightness among males and females and ventral and dorsal sides.

The file "dorso_ventral.R" runs the comparison of evolutionary rates between dorsal and ventral sides for UV coordinates and visible-light coordinates, and the PGLS between dorso-ventral distance in the UV and in visible light.

The file "sisterspecies_sympatry_divergence.R" runs the generalized linear model to test for divergence in UV coordinates for closely related sister species.

The file "specimen_age.R" runs the linear mixed model to test for the effect of specimen age on mean UV brightness of the wings.

The file "patternize_sis_sp.R" runs the patternize analysis on the sister species to compare ventral and dorsal divergence between species for males and females. However, to run this analysis you have to run the "generate_outlines.py" file first with python to produce same sized images with corresponding outlines in new folders within the images folder.

Within the data folder :
The color_sis_sp.txt contains putative color categorization for the sister species that are used just to set the final number of colours retained for the automatic recolorize colour segmentation and patternize analysis.
