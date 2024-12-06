# MAM_box_model

A box model driver for the Modal Aerosol Module (Liu et al., 2012; Liu et al., 2016). 


## Developers: 

```
   Richcard C. Easter (v1.0, original version)  
   Hui Wan            (v1.1-v1.2) 
   Jian Sun           (v1.1-v1.2) 
   Kai Zhang          (v1.1-v1.2) 
```

## Versions

### v1.0 - Original version by Richcard C. Easter

    - Initial development  
    - Same module/data structure as in E3SM
    - MAM code is the same as in E3SM version 1
    - Output in ascii text

### v1.1 - Improved I/O and namelist control 

    - Convergence test development
    - NetCDF output 
    - Namelist control to switch on/off individual processes 

### v1.2 - Current version 

    - Single box setup 
    - Namelist control for meteorological/initial condition 
    - Test case for wateruptake verification using bisection

### Version under development/testing

    - Use GCM netcdf output as input data 
    - Example data for selected ARM sites 
    - Precision control 
    - Test configurations 
    - Control for sulfur chemistry setup 

## Documentation 

To be updated.  

## Test cases

We have developed offline drivers for aerosol microphysical processes including water uptake, nucleation, condensation, coagulation, aging, and redistribution of particles among modes of different size ranges. In addition to performing unit tests for individual aerosol processes, the box model is also used to test numerical convergence for time integration of multiple processes and the impact of operator splitting. 

Please contact us for more details. 

## Contact

Kai Zhang (kai.zhang@pnnl.gov) 

## Reference 

- Liu, X., Easter, R. C., Ghan, S. J., Zaveri, R., Rasch, P., Shi, X., Lamarque, J.-F., Gettelman, A., Morrison, H., Vitt, F., Conley, A., Park, S., Neale, R., Hannay, C., Ekman, A. M. L., Hess, P., Mahowald, N., Collins, W., Iacono, M. J., Bretherton, C. S., Flanner, M. G., and Mitchell, D.: Toward a minimal representation of aerosols in climate models: description and evaluation in the Community Atmosphere Model CAM5, Geosci. Model Dev., 5, 709–739, https://doi.org/10.5194/gmd-5-709-2012, 2012. 
- Liu, X., Ma, P.-L., Wang, H., Tilmes, S., Singh, B., Easter, R. C., Ghan, S. J., and Rasch, P. J.: Description and evaluation of a new four-mode version of the Modal Aerosol Module (MAM4) within version 5.3 of the Community Atmosphere Model, Geosci. Model Dev., 9, 505–522, https://doi.org/10.5194/gmd-9-505-2016, 2016.  
- Easter, R. C., S. J. Ghan, Y. Zhang, R. D. Saylor, E. G. Chapman, N. S. Laulainen, H. Abdul-Razzak, L. R. Leung, X. Bian, and R. A. Zaveri (2004), MIRAGE: Model description and evaluation of aerosols and trace gases, J. Geophys. Res., 109, D20210, https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2004JD004571, doi:10.1029/2004JD004571.



