[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/camerianm/ExoEvo/master)


## Current status (2019/06/27):
* Extremely early version (development began 2019/04/19)
* Not yet recommended for use in research or educational settings

### Purpose and assumptions:
* Evolves mantle potential temperature, temperature-sensitive thermal parameters, and Rayleigh number over time
* Assumes an iron-free, magnesium-rich mantle; defaults to Mg-bearing endmembers of Mg/Fe bearing solid solution series.
* Does not incorporate melting processes that could change heat-transport efficiency and [convecting-layer thickness](https://doi.org/10.1089/ast.2017.1695)
* Assumes that Nusselt-Rayleigh scaling laws apply (more plausible in stagnant-lid regime)
* Assumes strictly temperature-dependent viscosity
* Estimates bulk thermal parameters from bulk mineralogy via weighted linear mixing/averaging of end-member properties
* User specifies basic planetary properties (# Earth masses, # Earth radii, Earth-normalized abundance of heat-producing elements, and bulk mineralogy); model estimates a plausible core mass fraction, core radius fraction, and bulk mantle density.
* Program interpolates specific heat capacity and thermal expansivity from P-T grids for each mineral end-member. Grids generated using [Stixrude and Lithgow-Bertelloni, 2011.](https://doi.org/10.1111/j.1365-246X.2010.04890.x)
* Direct importing of [ExoPlex](https://github.com/CaymanUnterborn/ExoPlex) / [PerpleX](http://www.perplex.ethz.ch/) outputs for broad-brush and radial-averaged properties (radius, core radius fraction, "representative" average pressure for thermal parameters, and bulk mantle mineralogy) 
* Allows easy addition of new minerals, by adding a dictionary entry (to mineralDB.py) with relevant entries.
* P/T grid for mineral's thermal parameters can be generated in ENKIportal via scripts provided in alphagrid_README.txt and CPgrid_README.txt. Users who want to add end-members, but who do not have ENKIportal access, should absolutely contact repository owner to request that a thermal grid be made.

### Not yet added:
   * Composition-specific thermal conductivity (in progress) **[~]**
   * Composition-specific temperature-dependent viscosity parameters (e.g. prefactors and activation energy for diffusion creep) **[~]**
   * Depth variations in T, P, and composition, to calculate evolution not of a 0-D mixture (bulk planet average at a single temperature and representative pressure, for each timestep), but a 1-D structure and accompanying adiabat

**[~]** = defaults to commonly-cited constant values for olivine
