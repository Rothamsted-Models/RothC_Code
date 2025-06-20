C******************************************************************************
C  Wrapper for RothC model 
C
C  February 2024
C
C  June 2025 this is the code that includes the Farina et al (2013) version of the model 
C
C  Farina et al, 2013, Geoderma. 200, 18-30, 10.1016/j.geoderma.2013.01.021
C  
C  Kevin Coleman
C
C******************************************************************************   
C
C INPUTS: 
C
C clay:        clay content of the soil (units: %)
C depth:       depth of topsoil (units: cm)
C IOM:         inert organic matter (t C /ha)
C nsteps:      number of timesteps 

C
C The following are needed for the Farina modification to the model (Farina et al, 2013, Geoderma. 200, 18-30, 10.1016/j.geoderma.2013.01.021)
C slit:        silt content of the soil (units: %) 
C BD:          bulk density (units: g/cm3)
C OC:          organic carbon (units: %)
C minRM_Moist: the minimum value the rate modifying factor for moisture can be (units: -, default=0.2)
C
C the following switches are needed to allow the user to choose which model option to run
C
C opt_RMmoist !  1: Standard RothC soil water parameters,
C             !  2: Van Genuchten soil properties and soil is allowed to be drier (ie hygroscopic / capillary water, -1000bar)
C             !  3: Van Genuchten soil properties, but uses the Standard RothC soil water function
C      
C opt_SMDbare !  1: Standard RothC bareSMD, 
C             !  2: bareSMD is set to wilting point -15bar (could be better for dry soils)
C
C year:     year
C month:    month (1-12)
C modern:   %modern 
C TMP:      Air temperature (C)
C Rain:     Rainfall (mm)
C Evap:     open pan evaporation (mm)
C C_inp:    Carbon input to the soil each month (units: t C /ha)
C FYM:      Farmyard manure input to the soil each month (units: t C /ha)
C PC:       Plant cover (0 = no cover, 1 = covered by a crop)
C DPM/RPM:  Ratio of DPM to RPM for carbon additions to the soil (units: none)

C  Note:
C  The shell reads in an example data set from RothC_input.dat, if your data is in another format you can change the read statements.
C
C  This model uses the first 12 months of weather (temp, rain, and evap), and land management information (C input, FYM input, and plant cover) to run to equilibrium    
C
      program RothC_shell
      
      implicit none
      
      integer MAXsteps

      parameter (MAXsteps = 73000)
      
      integer nsteps
      
      integer timeFact 
       
      integer i, j, k, k_month
      
      integer t_PC(MAXsteps)
      
      integer PC     
      
      integer t_year(MAXsteps)
      
      integer t_month(MAXsteps)
      
      real*8 t_mod(MAXsteps)
      
      real*8 t_tmp(MAXsteps)
      
      real*8 t_rain(MAXsteps)
      
      real*8 t_evap(MAXsteps)    
      
      real*8 t_C_Inp(MAXsteps)
      
      real*8 t_FYM_Inp(MAXsteps)
      
      real*8 t_DPM_RPM(MAXsteps)
      
      real*8 t_fert_N(MAXsteps)
      
      real*8 clay  ! clay content (Units: %)
      
      real*8 depth ! depth of topsoil (units: cm)
      
      real*8 silt  ! silt content (units: %) needed for the farina (2013) version
     
      real*8 BD    ! bulk density (units: g/cm3) needed for the farina (2013) version
               
      real*8 OC    ! organic carbon (units: %) needed for the farina (2013) version
      
      real*8 minRM_Moist ! (units: -, default=0.2) needed for the farina (2013) version
      
      
      real*8 DPM, RPM, BIO, HUM, IOM, SOC, total_CO2
      
      real*8 DPM_Rage,RPM_Rage,BIO_Rage,HUM_Rage,IOM_Rage,Total_Rage
      
      
      real*8 DPM_Delta, RPM_Delta, Bio_Delta, Hum_Delta, IOM_Delta
      real*8 Total_Delta
      
      integer YEAR, MONTH
      
      integer opt_RMmoist !  1: Standard RothC soil water parameters,
                          !  2: Van Genuchten soil properties and soil is allowed to be drier (ie hygroscopic / capillary water, -1000bar)
                          !  3: Van Genuchten soil properties, but uses the Standard RothC soil water function
      
      integer opt_SMDbare !  1: Standard RothC bareSMD, 
                          !  2: bareSMD is set to wilting point -15bar (could be better for dry soils)
      
      real*8 TEMP, RAIN, PEVAP
      
      real*8 RM_TMP, RM_Moist, RM_PC
      
      real*8 modernC
      
      ! C_inp (Carbon input to the soil each month units: t C /ha)
      
      real*8 C_inp, FYM_Inp, DPM_RPM  
      
      real*8 SMD      
      
      real*8 toc0, toc1, test
      
      real*8 time_begin, time_end
      
      call cpu_time (time_begin)
      
      IOM_Rage = 50000.0
      IOM_Delta = (exp(-50000/8035.0) - 1.0) * 1000.0
      
      
      DPM = 0.0
      RPM = 0.0
      BIO = 0.0
      HUM = 0.0
      IOM = 0.0
      
      total_CO2 = 0.0
      
      DPM_Rage = 0.0
      RPM_Rage = 0.0
      BIO_Rage = 0.0
      Hum_Rage = 0.0
      IOM_Rage = 50000.0  

C set initial soil water content (deficit) 
      SMD = 0.0
      
      minRM_Moist = 0.2  ! 0.2 is the default value for minRM_Moist for the farina (2013) version
C
C READ IN INPUT DATA: START
C
C read in RothC input data file: data will be passed from other programs at some point  
      open(11, file='RothC_input.dat', status='unknown')    
	read(11,*)               ! line is for info only 
	read(11,*)               ! line is for info only 
	read(11,*)               ! line is for info only 
	read(11,*)               ! line is for info only 
	read(11,*) opt_RMmoist, opt_SMDbare       
	read(11,*)               ! line is for info only  
	read(11,*)               ! line is for info only 
      if (opt_RMmoist.eq.1)then
        read(11,*)clay, depth, iom, nsteps
      else
        read(11,*)clay, depth, iom, nsteps, silt, BD, OC, minRM_Moist
      endif
      read(11,*)               ! line is for info only 
      read(11,*)               ! line is for info only 

	do i = 1, nsteps
	  read(11,*)t_year(i), t_month(i), t_mod(i), t_tmp(i),t_rain(i),
     &      t_evap(i), t_C_Inp(i), t_FYM_Inp(i), t_PC(i), t_DPM_RPM(i)
      enddo
      
      close(11)   
      
      
      open(71, file='year_results.out', status ='unknown')     
      open(91, file='month_results.out', status ='unknown')
C
C READ IN INPUT DATA: END
C      
      
C      
C run RothC to equilibrium: START
C uses first 12 months of weather and carbon data 
C
      k = 0
      j = 0
      
      
      SOC = DPM+RPM+Bio+Hum+IOM
      
        write(71,7100)
7100  format(5x, 'Year,', 2x,  'Month,', 
     &  1x, 'DPM_t_C_ha,', 1x,  'RPM_t_C_ha,', 
     &  1x, 'BIO_t_C_ha,', 1x, 'HUM_t_C_ha,',
     &  1x, 'IOM_t_C_ha,', 1x, 'SOC_t_C_ha,', 
     &  1x, 'CO2_t_C_ha,', 1x, ' deltaC')     
             
      write(71,101) j, DPM, RPM, Bio, Hum, iom, SOC,total_CO2
101   format(1x, '       0,', 1x, i6, ',', 7(f11.4, ','), ' -998.02')  

      timeFact = 12 ! monthly
          
      test = 100.0
      
       write(91,9100)
9100  format(4x, 'Year,',1x,  'Month,',1x, 'C_Inp_t_C_ha,', 
     &  1x,  'FYM_Inp_t_C_ha,', 1x,  'TEMP_C,', 1x, 'RM_TMP,',
     &  1x, 'RAIN_mm,', 1x, 'PEVAP_mm,',1x, 'SMD_mm,',
     &  1x,'RM_Moist,', 1x, 'PC,', 1x,  'RM_PC,',  
     &  1x, 'DPM_t_C_ha,', 1x,  'RPM_t_C_ha,', 
     &  1x, 'BIO_t_C_ha,', 1x, 'HUM_t_C_ha,',
     &  1x, 'IOM_t_C_ha,', 1x, 'SOC_t_C_ha,', 
     &  1x, 'CO2_t_C_ha') 
     
      YEAR = t_year(1)
      
      write(91,9101) Year, j, DPM, RPM, BIO, HUM, IOM, SOC, total_CO2
     
9101     format(1x, i7, ',', i6, ',', 13x ',', 15x, ',', 7x, ',',
     &         7x, ',', 8x, ',', 9x, ',', 7x, ',', 9x, ',',
     &           3x, ',', 6x, ',',f11.4, ',',f11.4,',', f11.4, ',',
     &        f11.4, ',',f11.4, ',',f11.4, ',',f11.4)   
         
      
      do ! Run to equililibrium: cycles through the first 12 months
       k = k + 1
       j = j + 1 
       
       if(k.eq.timeFact+1)k = 1   ! 13 if monthly
         if (test < 1E-6) exit
         YEAR = t_year(k)
         TEMP = t_tmp(k)
         RAIN = t_rain(k)
         PEVAP = t_evap(k)
         
         PC = t_PC(k)
         DPM_RPM = t_DPM_RPM(k)
         
         C_inp = t_C_Inp(k)
         FYM_Inp = t_FYM_Inp(k)
         
         modernC = t_mod(k) / 100.0             
         
         call RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC, total_CO2, 
     &     DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage, 
     &     modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM,
     &     C_Inp, FYM_Inp, SMD, RM_TMP, RM_Moist, RM_PC, 
     &     opt_RMmoist, opt_SMDbare, silt, BD, OC, minRM_Moist)  
        
         if(mod(k, timeFact)== 0)then 
           TOC0 = TOC1
           TOC1 = DPM+RPM+Bio+Hum
           test = abs(TOC1-TOC0)            
         endif    
         
      enddo
      
      total_CO2 = 0.0 ! reset CO2 to zero after the equilibrium run
      
      write(91,9102) Year, j-1, DPM, RPM, BIO, HUM, IOM, SOC, total_CO2
     
9102     format(1x, i7, ',', i6, ',', 13x ',', 15x, ',', 7x, ',',
     &         7x, ',', 8x, ',', 9x, ',', 7x, ',', 9x, ',',
     &           3x, ',', 6x, ',',f11.4, ',',f11.4,',', f11.4, ',',
     &        f11.4, ',',f11.4, ',',f11.4, ',',f11.4)   
     

C      
C run RothC to equilibrium: END
C
         
      Total_Delta = (exp(-Total_Rage/8035.0) - 1.0) * 1000.0   
      
      write(71,102) year, j-1, DPM, RPM, Bio, Hum, iom, SOC, total_CO2, 
     &              Total_Delta
102   format(1x, i8, ',', 1x, i6, ',', 7(f11.4,','),  f8.2)     
                                                             
C      
C run RothC for months 13 to the end: START
C     
      k_month = 0
      do i = timeFact+1, nsteps, 1   ! 13 if monthly
      
	  k_month = k_month + 1
            
        if(k_month.eq.timeFact+1)k_month = 1   ! 13 if monthly
         
        
         YEAR = t_year(i)
         TEMP = t_tmp(i)
         RAIN = t_rain(i)
         PEVAP = t_evap(i)
         
         PC = t_PC(i)
         DPM_RPM = t_DPM_RPM(i)
         
         C_inp = t_C_Inp(i)
         FYM_Inp = t_FYM_Inp(i)
         
         modernC = t_mod(i) / 100.0
           
         call RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC, total_CO2, 
     &     DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage, 
     &     modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM,
     &     C_Inp, FYM_Inp, SMD, RM_TMP, RM_Moist, RM_PC, 
     &     opt_RMmoist, opt_SMDbare, silt, BD, OC, minRM_Moist)    
         
         Total_Delta = (exp(-Total_Rage/8035.0) - 1.0) * 1000.0
         
         
         write(91,9103) Year, k_month, C_Inp, FYM_Inp, TEMP,RM_TMP, 
     &        RAIN, PEVAP, SMD, RM_Moist, PC, RM_PC,
     &        DPM,RPM,BIO,HUM, IOM, SOC, total_CO2
     
9103     format(1x, i7, ',', i6, ',', f13.3, ',',f15.3, ',',f7.1, ',',
     &         f7.4, ',',f8.1, ',',f9.1, ',',f7.2, ',', f9.4, ',',
     &           i3, ',',f6.1, ',',f11.4, ',',f11.4,',', f11.4, ',',
     &        f11.4, ',',f11.4, ',',f11.4, ',',f11.4)        


      if(mod(i, timeFact)== 0)then     ! print out results once a year
        write(71,103) year, DPM, RPM, Bio, Hum, IOM, SOC, total_CO2, 
     &                Total_Delta
      endif
         
      enddo  
C      
C run RothC for months 13 to the end: END
C   
      
      call cpu_time (time_end)
      
!      write(81,*) 'Time of operation was ', 
!     $    time_end - time_begin, ' seconds'
           
 103  format(1x,  i7, ',', 6x, '12,', 7(f11.4, ','), f8.2)  
 
      close (71)
      close (91)
  
      stop
      
      end
      
            