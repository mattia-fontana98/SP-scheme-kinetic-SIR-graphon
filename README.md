# README

**A. Bondesan, J. Borsotti, and M. Fontana**

This document is a guide to the MATLAB code used to reproduce the numerical results presented in the article "Kinetic models of opinion-driven epidemic dynamics modulated by graphons" co-authored by A. Bondesan, J. Borsotti, and M. Fontana (DOI: 10.48550/arXiv.2604.10614).

The repository is organized around a main script and a set of functions. The main script `kineticSIR_main.m` sets parameters, initializes data, runs the time loop, and saves outputs. All simulation settings are controlled from the first two sections of the script, namely "Model parameters" and "Control parameters".

The code was tested using MATLAB R2025b and requires no additional toolboxes.

## Model parameters

In this section, the user specifies the discretization in $(x, w)$, chooses the number of time steps `Nt` and the relaxation time `tau`, the drift and diffusion coefficients of the opinion dynamics `lambda` and `sigma_S`, `sigma_I`, `sigma_R`, and the epidemiological coefficients `beta_contact`, `gamma_recovery` together with the exponent `alpha`.

## Control parameters

In this section, the user selects the numerical options such as the quadrature rule for the Chang–Cooper weights through `scheme_w`, the time integration for the Fokker–Planck step via `scheme_t`, the initial condition `choice_IC`, and the functional forms of the graphon $\mathcal{B}$ and interaction functions $P$ and $G$ through `choice_B`, `choice_P`, and `choice_G`.

The flags `choice_B`, `choice_P`, and `choice_G` select among several predefined model variants. The corresponding parameters $\xi$, $r$, $\chi$, $a$ for $\mathcal{B}$ and $P$, and $\Delta$ for $G$ are kept inside the definition files; they can be modified directly in `compute_B.m`, `compute_P.m`, and `compute_G.m`.

At the end of each run, `kineticSIR_main.m` automatically saves the full workspace to `Output_kineticSIR.mat`, which is used by the plotting script. In addition, the flag `save_data` controls whether a separate `.mat` file containing only the arrays `x`, `w`, `ftot_over_time`, and `time` is written to disk. Before running the kinetic SIR simulation, the user must set `save_data = 1` if a subsequent popularity simulation is intended, as `popularity_main.m` loads this file as its input.

## Popularity module

Simulations for the popularity model are provided in the subfolder `Popularity_spread/`. This module is run separately from the kinetic SIR code and is driven by its own main script, `popularity_main.m`. All relevant settings are collected in the first part of the script, where the user can modify the model and the numerical parameters.

The popularity model is not standalone, as it requires as input the output of a previously computed kinetic SIR simulation. In particular, `popularity_main.m` loads the time-dependent distribution `ftot_over_time`, the corresponding time grid `time`, and the spatial–opinion grids `x`, `w` produced by `kineticSIR_main.m`.

In practice, the user should first enable the saving section in `kineticSIR_main.m`, so that the required arrays are written to disk. The popularity script then reads that saved dataset through a `load(...)` instruction. Therefore, the popularity dynamics is coupled to the previously computed kinetic SIR simulation and cannot be launched independently without such an input file.

Within `popularity_main.m`, the user can modify the corresponding model parameters, namely the popularity drift and diffusion coefficients `mu` and `zeta`, the spread factor `theta`, and the hypothesized danger level `w_hat`. The functional form of the propensity to interact $p$ is selected through `choice_p` and must be coherent with the choice made in `kineticSIR_main.m`, while the parameters $\xi$, $r$, and $\chi$ are kept inside the definition file; they can be modified in `compute_p.m`. The numerical options, similarly to those in `kineticSIR_main.m`, are controlled through `scheme_v` and `scheme_t`.

## Plotting

Both the kinetic SIR code and the popularity module are accompanied by dedicated plotting directories, `Plot_kineticSIR/` and `Popularity_spread/Plot_popularity/` respectively. Each folder contains a main plotting script, `plot_kineticSIR.m` and `plot_popularity.m` respectively, together with helper functions handling figure creation, style settings, and snapshot selection. The plotting scripts load the corresponding output files, namely `Output_kineticSIR.mat` and `Output_popularity.mat`, from the parent directory. The user can select which figures to generate through the switches at the top of each plotting script.
