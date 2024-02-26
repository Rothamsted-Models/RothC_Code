# The Rothamsted carbon model (RothC)

## Purpose

Roth C models the turnover of organic carbon in non-waterlogged top-soil.  It accounts for the effects of soil texture, temperature, moisture content and plant cover on the turnover process. It uses a monthly time step to calculate total organic carbon (t ha<sup>-1</sup>), microbial biomass carbon (t ha<sup>-1</sup>) and Δ<sup>14</sup>C (from which the equivalent radiocarbon age of the soil can be calculated). 

## Development history

The first version of RothC created by David Jenkinson and James Rayner in 1977 (Jenkinson and Rayner, 1977).

In 1987 an updated version was published, see Jenkinson et al. (1987).  This version included the prediction of the radiocarbon age of the soil, the pools POM (physically stabilized organic matter) and COM (chemically stabilized organic matter) were replaced with Hum (humified organic matter) and IOM (inert organic matter), and the microbial biomass pool was split into BioA (autochthonous biomass) and BioZ (zymogenous biomass).  

**In 1990, the two biomass pools were combined into a single pool (Jenkinson, 1990) this version is the standard version of the model, that this code refers to.**

Other published developments of the model include:

Farina et al. (2013) modified the soil water dynamics for semi-arid regions.

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
3.	created monthly and yearly outputs
The file can be used to create a standalone exe, or you can replace it with your own code to call and run RothC. Details of the inputs required, pools modelled, and units are in the code.


### RothC_input.dat  
This file contains input variables for the model.  

At the start of the file values for **clay** (%), **soil depth** (cm), **inert organic matter** (IOM, t C ha<sup>-1</sup>) and **number of steps** (nsteps) are recorded.  
Following that there is a table which records monthly data on **year**, **month**, **percentage of modern carbon**  (%), **mean air temperature** (Tmp, °C), **total monthly rainfall** (Rain, mm), **total monthly open-pan evaporation** (Evap, mm), **all carbon input entering the soil** (from plants, roots, root exudates) (C_inp, t C ha<sup>-1</sup>), **carbon input from farmyard manure** (FYM, t C ha<sup>-1</sup>), **plant cover** (PC, 0 for no plants e.g. bare or post-harvest, 1 for plants e.g. crop or grass), and the **DPM/RPM ratio** (DPM_RPM) of the carbon inputs from plants.

### year_results.out
This file contains the yearly values of the SOC (both the pools and Total) and the delta 14-carbon.

The pools are:  
**Year**  
**Month** 	- Always December for the yearly output  
**DPM** 	- Decomposable plant material (t C ha<sup>-1</sup>)  
**RPM** 	- Resistant plant material (t C ha<sup>-1</sup>)  
**BIO** 	- Microbial biomass (t C ha<sup>-1</sup>)  
**HUM**	- Humified organic matter (t C ha<sup>-1</sup>)  
**IOM** 	- Inert organic matter (t C ha<sup>-1</sup>)  
**SOC**	- Total soil organic carbon (t C ha<sup>-1</sup>)  
**deltaC** 	- delta <sup>14</sup>C (‰)  


The total organic carbon (soil organic carbon) is equal to the sum of the 5 pools. 

TOC or SOC = DRM + RPM + BIO + HUM + IOM 
     
### month_results.out
This file contains the monthly inputs, rate modifying factors, SOC pools.

**Year**  
**Month**  
**C_Inp_t_C_ha**		- C input (t C ha<sup>-1</sup>)  
**FYM_Inp_t_C_ha**	- Farmyard manure (t C ha<sup>-1</sup>)  
**TEMP_C**		- Air temperature (C)  
**RM_TMP**		- Rate modifying factor for temperature (-)  
**RAIN_mm**		- Rainfall (mm)  
**PEVAP_mm**		- Open pan evaporation (mm)  
**SWC_mm**		- Accumulated soil water deficit (mm)  
**RM_Moist**		- Rate modifying factor for soil moisture (-)  
**PC**			- Soil plant cover (0 bare or 1 covered)  
**RM_PC**			- rate modifying factor for crop cover  
**DPM_t_C_ha**		- Decomposable plant material (t C ha<sup>-1</sup>)  
**RPM_t_C_ha**		- Resistant plant material (t C ha<sup>-1</sup>)  
**BIO_t_C_ha**		- Microbial biomass (t C ha<sup>-1</sup>)  
**HUM_t_C_ha**		- Humified organic matter (t C ha<sup>-1</sup>)  
**IOM_t_C_ha**		- Inert organic matter (t C ha<sup>-1</sup>)  
**SOC_t_C_ha**		- Total soil organic carbon (t C ha<sup>-1</sup>)  

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

