[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/camerianm/ExoEvo/master)


## Current status (December 8, 2019):
* Extremely early version (development began April 2019)
* Not yet recommended for use in research or educational settings

### Purpose and assumptions:
* Evolves mantle potential temperature, temperature-sensitive thermal parameters, and Rayleigh number over time
* Assumes that Nusselt-Rayleigh scaling laws apply (more plausible in stagnant-lid regime)
* Assumes an iron-free, magnesium-rich mantle; defaults to Mg-bearing endmembers of Mg/Fe bearing solid solution series.
* Does not incorporate melting processes that could change heat-transport efficiency or [convecting-layer thickness](https://doi.org/10.1089/ast.2017.1695). 
* Assumes strictly temperature-dependent viscosity, diffusion creep (n=1) as deformation mode
* User specifies basic planetary properties (# Earth masses, # Earth radii, Earth-normalized abundance of heat-producing elements, and bulk mineralogy); model estimates a plausible core mass fraction, core radius fraction, and bulk mantle density.
* Without ExoPlex files, estimates bulk thermodynamic and physical parameters from bulk mineralogy via weighted linear mixing/averaging of end-member properties
* 'dynamic' mode interpolates specific heat capacity and thermal expansivity at each temperature the model traverses, using P-T grids for each mineral end-member. Grids generated using ENKIPortal / ThermoEngine implementation of [Stixrude and Lithgow-Bertelloni, 2011.](https://doi.org/10.1111/j.1365-246X.2010.04890.x)
* Adaptive to translation issues between deformation data and numerical settings - i.e., the model can work with **either** scaled prefactors (i.e. a reference viscosity 'visc0' in Pa s, at a reference temperature 'scaletemp' in K) and unscaled prefactors for diffusion creep. NOTE: this currently operates through a hard-coded threshhold, where **if planet['visc0']>1.0e13,** planet['visc0'] is assumed to be reference viscosity at a scaling temperature planet['scaletemp'].
* Permits direct importing of [ExoPlex](https://github.com/CaymanUnterborn/ExoPlex) / [PerpleX](http://www.perplex.ethz.ch/) CSV-format outputs for (a) structural parameters and (b) broad-brush, volume-averaged or mass-averaged properties (alpha, Cp, k, and mass percent of mineral species).
* Allows easy addition of new minerals, by adding a dictionary entry (to mineralDB.py) with relevant entries, as well as P/T grid of properties from ENKIPortal
* P/T grid for mineral's thermal parameters can be generated in ENKIportal via scripts provided in alphagrid_README.txt and CPgrid_README.txt. Users who want to add end-members, but who do not have ENKIportal access, should contact repository owner to request that a thermal grid be made.

### Not yet added:
   * Composition-specific thermal conductivity (in progress) **[~]**
   * Composition-specific temperature-dependent viscosity parameters (e.g. prefactors and activation energy for diffusion creep) **[~]**
   * Depth variations in T, P, and composition, to calculate evolution not of a 0-D mixture (bulk planet average at a single temperature and representative pressure, for each timestep), but a 1-D structure and accompanying adiabat

**[~]** = defaults to commonly-cited constant values for olivine
