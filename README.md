### Replication Code for Paper on an Incentive-Compatible Draft Allocation Mechanism

This repository contains Julia and R code for replicating the figures and tables in 

The `data` folder contains cleaned data for the NBA from 1985 to 1989, and from
2015 to 2019.

The code for replicating the figures and tables in the papers is organized as follows:

- `draftpolicy.jl` contains functions for calculating R-NTD on real and simulated data
- `nba_sims.jl` contains code to simulate the performance of R-NTD
