# Documentation

[![N|Solid](https://cldup.com/dTxpPi9lDf.thumb.png)](https://nodesource.com/products/nsolid)

This repository contains the source code for the following paper, please pursue [this link](http://www.sciencedirect.com/science/article/pii/S0360544217310423?via%3Dihub):

Li, B., de Queiroz, A. R., DeCarolis, J. F., Bane, J., He, R., Keeler, A. G., Neary, V. S. The economics of electricity generation from Gulf Stream currents. Energy, 134 (2017): 649 - 658.

# Function definitions
- `MEP_calculate.m`
Input : `<year>.mat`, `depth_domain.mat`
Output: Save monthly energy output (`MEP_grid`) in `MEP<year>.mat`, each year is an independent file

- `MEP_combine.m`
Input: `MEP<year>.mat`
Output: A single file (`MEP.mat`) including 6 years of monthly energy output (`MEP`)

- `pa_preprocess.m`
Input: `MEP.mat`, `2009result.mat`
Output: [Figure 3](http://www.sciencedirect.com/science/article/pii/S0360544217310423?via%3Dihub#fig3) (coefficient of correlation vs. distance) and [Figure A](https://ars.els-cdn.com/content/image/1-s2.0-S0360544217310423-mmc1.pdf) (in the Appendix, frequency vs. distance).

- `pa.m`
This script runs the portfolio optimization as a continuous quadratic programming (QP) problem.
Input: `MEP.mat`, `depth_domain.mat`, `2009result.mat`
Output: A figure showing CF vs. 𝜎2.

- `pa_miqp.m`
This script runs the portfolio optimization as a mixed integer quadratic programing (MIQP) problem.
Input: `MEP.mat`, `depth_domain.mat`, `2009result.mat`, `<year>result.mat` if run 6-site test case
Output: Due to computational complexity, this model is never solved and a two-stage simplified version is solved instead.

- `pa_2step.m`
This script runs the portfolio optimization as a two-stage problem, the first stage is a continuous linear programming problem (LP) and the second step is a mixed integer quadratic programing (MIQP) problem.
Input: `MEP.mat`, `depth_domain.mat`, `2009result.mat`
Output: A scatter graph showing Pareto frontier of both stages
