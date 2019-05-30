[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/camerianm/ExoEvo/master)


## Current status (2019/05/30):
* Extremely early version (development began 2019/04/19)
* Not yet recommended for use in research or educational settings

### Purpose and assumptions:
* Evolves mantle potential temperature, temperature-sensitive thermal parameters, and Rayleigh number over time
* Assumes an iron-free, magnesium-rich mantle; defaults to Mg-bearing endmembers of Mg/Fe bearing solid solution series.
* Addition of new minerals is as simple as adding a dictionary entry (to mineralDB.py) with relevant coefficients to approximate curve from Berman, et al. (*NOTE: This process will be further simplified before July 2019.*)
* Does not incorporate melting processes that could change heat-transport efficiency and [convecting-layer thickness](https://doi.org/10.1089/ast.2017.1695)
* Assumes that Nusselt-Rayleigh scaling laws apply (more plausible in stagnant-lid regime)
* Assumes strictly temperature-dependent viscosity
* User specifies basic planetary properties (# Earth masses, # Earth radii, Earth-normalized abundance of heat-producing elements, and bulk mineralogy); model estimates a plausible core mass fraction, core radius fraction, and bulk mantle density.
* Estimates bulk thermal parameters from bulk mineralogy via weighted linear mixing/averaging
* Specific heat capacity and thermal expansivity parameters at 75 GPa from [Stixrude and Lithgow-Bertelloni, 2011.](https://doi.org/10.1111/j.1365-246X.2010.04890.x)
* Direct importing and radial-averaging functions for [ExoPlex](https://github.com/CaymanUnterborn/ExoPlex) / [PerpleX](http://www.perplex.ethz.ch/) output (composition only)

### Not yet added:
   * ExoPlex depth, pressure, and temperature radial profile import functions
   * P-T grid of thermal expansivity and Cp for each end-member mineral included
   * Composition-specific thermal conductivity **[~]**
   * Composition-specific temperature-dependent viscosity parameters (e.g. prefactors and activation energy for diffusion creep) **[~]**
   * Radial variations in T, P, and composition, to calculate evolution not of a 0-D mixture (bulk planet average at a single temperature for each timestep), but a 1-D structure and accompanying adiabat

**[~]** = defaults to commonly-cited constant values for olivine
