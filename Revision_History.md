# Revision History

## Version 2.0.0
This version adds options to align with the Farina et al. (2013) which explored changing the effect of soil moisture for semi-arid soils.
The model can be run as before (version 1.0.0) through specific options.

Options are for the soil moisture rate modifying factor: opt_RMmoist and opt_SMDbare.

When opt_RMmoist is given the value 2 or 3, four additional variables are required:
Silt (%), bulk density (cm<sup>-3</sup>), organic carbon (%), and the value to use as the minimum for soil moisture rate modifying factor (default value aligning with version 1.0.0 is 0.2).
These are not required when opt_RMmoist is given the value 1 (the Fortran code will not read them if present).

Results files now include accumulated CO<sub>2</sub> as a column.

The example files have been updated to show the new structures of the files.
See README.md and RothC_description.pdf for more details.

## Version 1.0.0
This version was the initial release of RothC aligning with Jenkinson (1990) (also previously referred to as RothC 26.3).