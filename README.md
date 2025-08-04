# The Rothamsted carbon model (RothC)

## Purpose

Roth C models the turnover of organic carbon in non-waterlogged top-soil.  It accounts for the effects of soil texture, temperature, moisture content and plant cover on the turnover process. It uses a monthly time step to calculate total organic carbon (t ha<sup>-1</sup>), microbial biomass carbon (t ha<sup>-1</sup>) and Δ<sup>14</sup>C (from which the equivalent radiocarbon age of the soil can be calculated). 

## Development history

The first version of RothC was created by David Jenkinson and James Rayner in 1977 (Jenkinson and Rayner, 1977).

In 1987 an updated version was published, see Jenkinson et al. (1987).  This version included the prediction of the radiocarbon age of the soil, the pools POM (physically stabilized organic matter) and COM (chemically stabilized organic matter) were replaced with Hum (humified organic matter) and IOM (inert organic matter), and the microbial biomass pool was split into BioA (autochthonous biomass) and BioZ (zymogenous biomass).  

**In 1990, the two biomass pools were combined into a single pool (Jenkinson, 1990) this version is the standard version of the model**
**Farina et al. (2013) modified the soil water dynamics for semi-arid regions.**
**The code from v2.0.0 includes the functionality from both of these developments**

Other published developments of the model include:

Giongo et al. (2020) created a daily version and modified the soil water dynamics, for Caatinga shrublands, in the semiarid region, North-East Brazil.

 
## Description of files included

### RothC_description.docx
This file contains the description of the model.


### RothC.for
This file contains the RothC code, it can be used as a standalone subroutine or used with shell.for to create an exe file. Details of the inputs required, pools modelled, and units are in the code.


### Shell.for
This file is intended as an example of how to: 
1.	read in the input data
2.	call the subroutine
3.	create monthly and yearly outputs
The file can be used to create a standalone exe, or you can replace it with your own code to call and run RothC. Details of the inputs required, pools modelled, and units are in the code.


### RothC_input.dat  
This file contains input variables for the model.  

At the start of the file are two options relating to the soil moisture function.
**If 1 is provided as the value for both then RothC will run as standard (like v1.0.0).**
**opt_RMmoist**: If 2 is provided then the soil can dry further than standard (i.e. -1500 kPa to -100000 kPa), without further reduction in RM_Moist
**opt_RMmoist**: If 3 is provided then the soil moisture function follows the same structure as standard but the moisture deficits corresponding to -100 kPa and -1500 kPa are calculated by Mualem-van Genuchten equations (van Genuchten, 1980; Wösten et al., 1999) and an adjustment is made to how the bare soil modification is calculated.
**opt_SMDbare**: If 2 is provided then the bare soil adjustment in the moisture function moves to the moisture deficit corresponding to -1500 kPa.

Next, values for **clay** (%), **soil depth** (cm), **inert organic matter** (IOM, t C ha<sup>-1</sup>) and **number of steps** (nsteps) are provided.
If 2 is provided for **Opt_RMmoist** then a further four variables are expected: **silt** (%), **bulk density** (g cm<sup>-3</sup>), **organic carbon** (%), and **minRM_Moist** which is the minimum value for the rate modifying factor for moisture (0.2 in standard version, and tested with 0.15 and 0.1 in (Farina et al. 2013))
**opt_RMmoist**: If 1 is provided then the additional four variables can still be present in the input file but will not be read in. 

Following those there is a table which provides monthly data of **year**, **month**, **percentage of modern carbon**  (%), **mean air temperature** (Tmp, °C), **total monthly rainfall** (Rain, mm), **total monthly open-pan evaporation** (Evap, mm), **all carbon input entering the soil** (from plants, roots, root exudates) (C_inp, t C ha<sup>-1</sup>), **carbon input from organic amendment** (OA, t C ha<sup>-1</sup>), **plant cover** (PC, 0 for no plants e.g. bare or post-harvest, 1 for plants e.g. crop or grass), the **allocations of plant material to DPM and RPM pools** (PL_DPM_f and PL_RPM_f; sum = 1), and **allocations of organic amendment to DPM, RPM, Bio, and Hum pools** (OA_DPM_f, OA_RPM_f, OA_Bio_f, and OA_Hum_f; sum = 1).

### year_results.out
This file contains the yearly values of the SOC (both the pools and Total) and the delta 14-carbon.

The pools are:  
**Year**  
**Month** 	    - Always December for the yearly output  
**DPM_t_C_ha** 	- Decomposable plant material (t C ha<sup>-1</sup>)  
**RPM_t_C_ha** 	- Resistant plant material (t C ha<sup>-1</sup>)  
**Bio_t_C_ha** 	- Microbial biomass (t C ha<sup>-1</sup>)  
**Hum_t_C_ha**	- Humified organic matter (t C ha<sup>-1</sup>)  
**IOM_t_C_ha** 	- Inert organic matter (t C ha<sup>-1</sup>)  
**SOC_t_C_ha**	- Total soil organic carbon (t C ha<sup>-1</sup>)
**CO2_t_C_ha**  - Accumulated CO<sub>2</sub> (t C ha<sup>-1</sup>)  
**deltaC** 	    - delta <sup>14</sup>C (‰)  


The total organic carbon (soil organic carbon) is equal to the sum of the 5 pools. 

TOC or SOC = DRM + RPM + Bio + Hum + IOM 

CO2 is set to 0 after the spin-up period and accumulates monthly as the sum of DPM_co2, RPM_co2, Bio_co2, and Hum_co2 each month.

### month_results.out
This file contains the monthly inputs, rate modifying factors, SOC pools.

**Year**  
**Month**  
**C_Inp_t_C_ha**    - C input (t C ha<sup>-1</sup>)  
**OA_Inp_t_C_ha**	- Farmyard manure (t C ha<sup>-1</sup>)  
**TEMP_C**		    - Air temperature (C)  
**RM_TMP**		    - Rate modifying factor for temperature (-)  
**RAIN_mm**		    - Rainfall (mm)  
**PEVAP_mm**		- Open pan evaporation (mm)  
**SWC_mm**		    - Accumulated soil water deficit (mm)  
**RM_Moist**		- Rate modifying factor for soil moisture (-)  
**PC**			    - Soil plant cover (0 bare or 1 covered)  
**RM_PC**			- rate modifying factor for crop cover  
**DPM_t_C_ha**		- Decomposable plant material (t C ha<sup>-1</sup>)  
**RPM_t_C_ha**		- Resistant plant material (t C ha<sup>-1</sup>)  
**Bio_t_C_ha**		- Microbial biomass (t C ha<sup>-1</sup>)  
**Hum_t_C_ha**		- Humified organic matter (t C ha<sup>-1</sup>)  
**IOM_t_C_ha**		- Inert organic matter (t C ha<sup>-1</sup>)  
**SOC_t_C_ha**		- Total soil organic carbon (t C ha<sup>-1</sup>)  
**CO2_t_C_ha**      - Accumulated CO<sub>2</sub> (t C ha<sup>-1</sup>) 

## Requirements
The code does not require any particular of version of Fortran, so can be compiled in both windows and Linux.

## Installation/set-up
A Fortran compiler is needed.


The code can be used in the following ways:
1.	The two Fortran files (RothC.for and Shell.for) can either be compiled and linked to create a standalone exe, which uses the input file (RothC_input.dat), when run, monthly (month_results.out) and yearly (year_results.out) output files are created.  
2.	The file (shell.for) can be modified to read in required data in the format you have, your modified code can be compiled and linked to RothC.for. 
3.	The file (RothC.for) can be called by your exiting code as a subroutine.    


**Example of how to run the model**  
The file RothC_input.dat contains all the inputs data needed to run the model. The month results (month_results.out) and year results (year_results.out) files correspond to this input file as an example. 
The model is normally run to equilibrium using average temperature, rainfall, open pan evaporation, an average carbon input to the soil, the equilibrium run is to initialise the soil carbon pools. Once the soil carbon pools have been initialised, the model is run for the period of interest. The met data (temperature, rainfall and evaporation) can be average or actual weather data. The carbon input to the soil can be: 1) adjusted so the modelled output matches the measured data, or 2) can be estimated from yield data (Bolinder et al., 2007), or NPP data.  


## References

Bolinder MA, Janzen HH, Gregorich EG, Angers DA, VandenBygaart AJ. An approach for estimating net primary productivity and annual carbon inputs to soil for common agricultural crops in Canada. Agriculture, Ecosystems & Environment 2007; 118: 29-42.  
Farina R, Coleman K, Whitmore AP. Modification of the RothC model for simulations of soil organic C dynamics in dryland regions. Geoderma 2013; 200: 18-30.  
Giongo V, Coleman K, Santana MD, Salviano AM, Olszveski N, Silva DJ, et al. Optimizing multifunctional agroecosystems in irrigated dryland agriculture to restore soil carbon - Experiments and modelling. Science of the Total Environment 2020; 725.  
Jenkinson DS. The Turnover of Organic-Carbon and Nitrogen in Soil. Philosophical Transactions of the Royal Society of London, Series B: Biological Sciences 1990; 329: 361-368.  
Jenkinson DS, Hart PBS, Rayner JH, Parry LC. Modelling the turnover of organic matter in long-term experiments at Rothamsted. INTECOL Bulletin 1987; 15: 1-8.  
Jenkinson DS, Rayner JH. Turnover of soil organic matter in some of the Rothamsted classical experiments. Soil Science 1977; 123: 298-305.  
van Genuchten, M.T., 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. Soil Science Society of America Journal 44 (5), 892–898.
Wösten, J.H.M., Lilly, A., Nemes, A., Le Bas, C., 1999. Development and use of a database of hydraulic properties of European soils. Geoderma 90 (3–4), 169–185.
