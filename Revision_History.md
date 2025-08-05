# Revision History

## Version 2.1.0
This version changes how the model receives the allocation of plant material and organic amendments to the active pools.
The proportions of plant material going to DPM and RPM (PL_DPM_f and PL_RPM_f, respectively; sum = 1) are specified monthly.
The proportions of organic amendment going to DPM, RPM, Bio, and Hum pools (OA_DPM_f, OA_RPM_f, OA_Bio_f, OA_Hum_f, respectively; sum = 1) are specified monthly.
In this way the plant material and organic amendment quality/composition can be specified more readily.

Previously, the plant material was split between DPM and RPM pools defined by a single value ratio (e.g. 1.44) provided in the input file and organic amendments were split between DPM, RPM, and Hum pools with values set within the model code.

Throughout the code: FYM (farmyard manure) has been replaced by OA (organic amendment); the case of the BIO or Bio and HUM or Hum pools have been made uniform as Bio and Hum (as these are more abbreviations than acronyms).

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