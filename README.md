# pvw

Simulation of the flux of nutrients through a viscous mucus layer to a spherical phytoplankton cell.  Mucus reduces both the molecular diffusivity of nutrients and the thermal conductivity.  A heat flux from the cell can be prescribed.  Both steady-state (pvws.m) and time-dependent (pvwt.m) codes are included.

## Requirements

These simulations require Matlab.  While open-source codebases are preferred, these were written at a time when my peers all used Matlab. The functions pvws.m and pvwt.m include the sub-functions SW_Conductivity and SW_Viscosity from the [Thermophysical Properties of Seawater](http://web.mit.edu/seawater) library v 3.1.2 (Nayar et al. 2016; Sharqawy et al. 2010).

## Usage

Create a `pvw` object using the pvw.m class.
```
par=pvw
```
This sets default values for all parameters used in pvws.m and pvwt.m. All parameters and their units are described in pvw.m lines 4-15.  As an example, adjust the parameters as follows:
```
par.scale=10;   % Set the cell size to 10 microns
par.b=0.1;      % Set the mucus layer thickness to ~10 microns
par.Vf=2;       % Set the mucus viscosity to 2 x seawater
```

For a steady-state simulation, use pvws.m:
```
pvws(par)
```
This assigns output variables to the base workspace.

To simulate a time-dependent nutrient pulse, use pvwt.m.  Additional pulse parameters in pvw.m lines 18-26 should be adjusted.
```
pvwt(par)
```
This saves output variables to a file with a unique name based on the run name `par.rn` and the current date and time.

## References

Nayar, K. G., Sharqawy, M. H., & Banchik, L. D. 2016. Thermophysical properties of seawater: A review and new correlations that include pressure dependence. Desalination, 390, 1-24.

Sharqawy, M. H., Lienhard, J. H., & Zubair, S. M. 2010. Thermophysical properties of seawater: a review of existing correlations and data. Desalination and water Treatment, 16(1-3), 354-380.
