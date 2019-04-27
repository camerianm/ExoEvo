[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/camerianm/ExoEvo/master?urlpath=https%3A%2F%2Fgithub.com%2Fcamerianm%2FExoEvo%2Fblob%2Fmaster%2FExoEvo.ipynb)

## Current status (2019/04/25):
* Extremely early version (development began 2019/04/19)
* Not yet recommended for use in research or educational settings

### Purpose and assumptions:
* Evolves mantle potential temperature, temperature-sensitive thermal parameters, and Rayleigh number over time
* Does not incorporate melting processes that could change heat-transport efficiency and [convecting-layer thickness](https://doi.org/10.1089/ast.2017.1695)
* Assumes that Nusselt-Rayleigh scaling laws apply (more plausible in stagnant-lid regime)
* Assumes strictly temperature-dependent viscosity
* User specifies basic planetary properties (# earth masses, # earth radii, core mass fraction, core radius fraction, ratio of heat-producing elements, and dominant mantle mineral)
* Specific heat capacity parameters from [Berman, et al.](https://doi.org/10.4095/223425) and citations therein

### Not yet added:
   * Composition-specific thermal conductivity **[~]**
   * Composition-specific thermal expansivity **[~]**
   * Composition-specific temperature-dependent viscosity parameters (e.g. prefactors and activation energy for diffusion creep) **[~]**
   * Weighted linear mixing of mineral phases (e.g. 70% Fo, 30% Fa)
   * Direct import function for [ExoPlex](https://github.com/CaymanUnterborn/ExoPlex) / [PerpleX](http://www.perplex.ethz.ch/) output
   * Predominant lower-mantle minerals (perovskite, bridgmanite, wustite) and transition-zone phases (majorite, ringwoodite, wadsleyite)

**[~]** = defaults to commonly-cited constant values for olivine


### Current mineral options:
* forsterite (default)
* fayalite
* orthoenstatite
* clinoenstatite
* periclase
* corundum
* spinel
* diopside
* diamond
* ca-al pyroxene