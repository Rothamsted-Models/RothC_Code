C******************************************************************************
C  Wrapper for RothC model 
C
C  October 2025
C
C  Kevin Coleman
C  Jonah Prout
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
C silt:        silt content of the soil (units: %) 
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
C opt_tstep !  1: Monthly time step
C           !  2: Daily time step
C
C opt_spin  !  1: use spin up
            !  2: initialize but read in average weather
            !  3: initialize but don't read in average weather
C
C year:     year
C tstep:    month (1-12) or day (1 to 365) depending on opt_tstep
C modern:   %modern 
C TMP:      Air temperature (C)
C Rain:     Rainfall (mm)
C Evap:     open pan evaporation (mm)
C Pl_inp:   Carbon input to the soil each month from plants (units: t C /ha)
C OA_inp:   Organic amendment input to the soil each month (units: t C /ha)
C PC:       Plant cover (0 = no cover, 1 = covered by a crop)
C Pl_DPM_f: Fraction of plant carbon to DPM
C Pl_DPM_f: Fraction of plant carbon to RPM
C OA_DPM_f: Fraction of organic amendment carbon to DPM
C OA_RPM_f: Fraction of organic amendment carbon to RPM
C OA_Bio_f: Fraction of organic amendment carbon to Bio
C OA_Hum_f: Fraction of organic amendment carbon to Hum

C  Note:
C  The shell reads in an example data set from RothC_input.dat, if your data is in another format you can change the read statements.
C
C  This model uses the first 12 months of weather (temp, rain, and evap), and land management information (C input, OA input, and plant cover) to run to equilibrium    
C
      program RothC_shell
      
      implicit none
      
      integer MAXsteps

      parameter (MAXsteps = 73000)
      
      integer nsteps
      
      integer i, j, k, k_tstep
      
      integer t_PC(MAXsteps)
      
      integer PC     
      
      integer t_year(MAXsteps)
      
      integer t_tstep(MAXsteps)
      
      real*8 t_mod(MAXsteps)
      
      real*8 t_tmp(MAXsteps)
      
      real*8 t_rain(MAXsteps)
      
      real*8 t_evap(MAXsteps)    
      
      real*8 t_Pl_inp(MAXsteps)
      
      real*8 t_OA_inp(MAXsteps)
      
      real*8 t_Pl_DPM_f(MAXsteps) ! fraction of plant C to DPM
      
      real*8 t_Pl_RPM_f(MAXsteps) ! fraction of plant C to RPM
      
      real*8 t_OA_DPM_f(MAXsteps) ! fraction of org amd C to DPM
      
      real*8 t_OA_RPM_f(MAXsteps) ! fraction of org amd C to RPM
      
      real*8 t_OA_Bio_f(MAXsteps) ! fraction of org amd C to Bio
      
      real*8 t_OA_Hum_f(MAXsteps) ! fraction of org amd C to Hum
      
      
      real*8 clay  ! clay content (Units: %)
      
      real*8 depth ! depth of topsoil (units: cm)
      
      real*8 silt  ! silt content (units: %) needed for the farina (2013) version
     
      real*8 BD    ! bulk density (units: g/cm3) needed for the farina (2013) version
               
      real*8 OC    ! organic carbon (units: %) needed for the farina (2013) version
      
      real*8 minRM_Moist ! (units: -, default=0.2) needed for the farina (2013) version
      
      real*8 DPM_init, RPM_init, Bio_init, Hum_init ! initial values of dpm, rpm, bio, and hum if not using spin up mode
      
      real*8 DPM, RPM, Bio, Hum, IOM, SOC, total_CO2
      
      real*8 DPM_Rage,RPM_Rage,Bio_Rage,Hum_Rage,IOM_Rage,Total_Rage
      
      
      real*8 DPM_Delta, RPM_Delta, Bio_Delta, Hum_Delta, IOM_Delta
      real*8 Total_Delta
      
      integer YEAR, TIMESTEP
      

      
      integer opt_RMmoist !  1: Standard RothC soil water parameters,
                          !  2: Van Genuchten soil properties and soil is allowed to be drier (ie hygroscopic / capillary water, -1000bar)
                          !  3: Van Genuchten soil properties, but uses the Standard RothC soil water function
      
      integer opt_SMDbare !  1: Standard RothC bareSMD, 
                          !  2: bareSMD is set to wilting point -15bar (could be better for dry soils)
      
      integer opt_tstep !  1: Monthly, tstep = 1/12
                        !  2: Daily, tstep = 1/365
      
      integer opt_Spin  !  1: use spin up
                        !  2: initialize but read in average weather
                        !  3: initialize but don't read in average weather
      
      integer start_loop      
      integer year_end
      
      real*8 TEMP, RAIN, PEVAP
      
      real*8 RM_TMP, RM_Moist, RM_PC
      
      real*8 modernC
      
      ! C_inp (Carbon input to the soil each tstep units: t C /ha)
      
      real*8 Pl_inp, OA_inp
      
      ! pool_f are the decimal fractions of plant and amendment C
      real*8 Pl_DPM_f, Pl_RPM_f, OA_DPM_f, OA_RPM_f, OA_Bio_f, OA_Hum_f
      
      real*8 SMD      
      
      real*8 toc0, toc1, test
      
      real*8 time_begin, time_end
      
      call cpu_time (time_begin)
      
      DPM = 0.0
      RPM = 0.0
      Bio = 0.0
      Hum = 0.0
      IOM = 0.0
      
      total_CO2 = 0.0
      
      DPM_Rage = 0.0
      RPM_Rage = 0.0
      Bio_Rage = 0.0
      Hum_Rage = 0.0
      IOM_Rage = 50000.0 
      
      IOM_Delta = (exp(-50000/8035.0) - 1.0) * 1000.0      

C set initial soil water content ( soil moisture deficit) 
      SMD = 0.0
      
      minRM_Moist = 0.2  ! 0.2 is the default value for minRM_Moist for the farina (2013) version
C
C READ IN INPUT DATA: START
C
C NOTE: Input file is different if opt_spin = 3 is choosen
C
C read in RothC input data file: data will be passed from other programs at some point  
      open(11, file='RothC_input.dat', status='unknown')    
	read(11,*)               ! line is for info only 
	read(11,*)               ! line is for info only 
	read(11,*)               ! line is for info only 
	read(11,*)               ! line is for info only 
	read(11,*) opt_RMmoist, opt_SMDbare, opt_tstep, opt_spin 
	read(11,*)               ! line is for info only  
	read(11,*)               ! line is for info only 
      if(opt_spin.eq.1)then
        read(11,*) iom
      else 
        read(11,*) iom, dpm_init, rpm_init, Bio_init, Hum_init
      endif 
      read(11,*)               ! line is for info only  
	read(11,*)               ! line is for info only 
      if (opt_RMmoist.eq.1)then
        read(11,*)nsteps, clay, depth
      else
        read(11,*)nsteps, clay, depth, silt, BD, OC, minRM_Moist
      endif
      read(11,*)               ! line is for info only 
      read(11,*)               ! line is for info only 

	do i = 1, nsteps
	  read(11,*)t_year(i), t_tstep(i), t_mod(i), t_tmp(i),t_rain(i),
     &      t_evap(i), t_Pl_inp(i), t_OA_inp(i), t_PC(i),
     &      t_Pl_DPM_f(i), t_Pl_RPM_f(i), t_OA_DPM_f(i), t_OA_RPM_f(i),
     &      t_OA_Bio_f(i), t_OA_Hum_f(i)
      enddo
      
      close(11)   
      
      
      open(71, file='year_results.out', status ='unknown')     
      open(91, file='tstep_results.out', status ='unknown')
C
C READ IN INPUT DATA: END
C      
      
      if(opt_tstep == 1)then
          year_end = 12
      elseif(opt_tstep == 2)then
          year_end = 365
      endif
      
C      
C run RothC to equilibrium: START
C uses first year of weather and carbon data 
C
      k = 0
      j = 0
      
      
      SOC = DPM+RPM+Bio+Hum+IOM
      
      if(opt_spin==1)then
        write(71,7100)
      else
        write(71,7101)
      endif
      
7100  format(5x, 'Year,', 2x,  'tstep,', 
     &  1x, 'DPM_t_C_ha,', 1x,  'RPM_t_C_ha,', 
     &  1x, 'Bio_t_C_ha,', 1x, 'Hum_t_C_ha,',
     &  1x, 'IOM_t_C_ha,', 1x, 'SOC_t_C_ha,', 
     &  1x, 'CO2_t_C_ha,', 1x, ' deltaC')     
             
7101  format(5x, 'Year,', 2x,  'tstep,', 
     &  1x, 'DPM_t_C_ha,', 1x,  'RPM_t_C_ha,', 
     &  1x, 'Bio_t_C_ha,', 1x, 'Hum_t_C_ha,',
     &  1x, 'IOM_t_C_ha,', 1x, 'SOC_t_C_ha,', 
     &  1x, 'CO2_t_C_ha,')   

      if(opt_spin==1)then
        write(71,101) j, DPM, RPM, Bio, Hum, iom, SOC,total_CO2
      else
        write(71,111) j, DPM, RPM, Bio, Hum, iom, SOC,total_CO2
      endif     
             
101   format(1x, '       0,', 1x, i6, ',', 7(f11.4, ','), ' -998.02')  
111   format(1x, '       0,', 1x, i6, ',', 7(f11.4, ','))  

      
       write(91,9100)
9100   format(4x, 'Year,',1x,  'tstep,',1x, 'Pl_inp_t_C_ha,', 
     &  1x,  'OA_inp_t_C_ha,', 1x,  'TEMP_C,', 1x, 'RM_TMP,',
     &  1x, 'RAIN_mm,', 1x, 'PEVAP_mm,',1x, 'SMD_mm,',
     &  1x,'RM_Moist,', 1x, 'PC,', 1x,  'RM_PC,',  
     &  1x, 'DPM_t_C_ha,', 1x,  'RPM_t_C_ha,', 
     &  1x, 'Bio_t_C_ha,', 1x, 'Hum_t_C_ha,',
     &  1x, 'IOM_t_C_ha,', 1x, 'SOC_t_C_ha,', 
     &  1x, 'CO2_t_C_ha') 
       
      test = 100.0   
      
      if(opt_spin==3)then
        year = 1
      else
        YEAR = t_year(1)
      endif
      
      write(91,9101) Year, j, DPM, RPM, Bio, Hum, IOM, SOC, total_CO2
     
9101     format(1x, i7, ',', i6, ',', 13x, ',', 15x, ',', 7x, ',',
     &         7x, ',', 8x, ',', 9x, ',', 7x, ',', 9x, ',',
     &           3x, ',', 6x, ',',f11.4, ',',f11.4,',', f11.4, ',',
     &        f11.4, ',',f11.4, ',',f11.4, ',',f11.4)   
         

      if (opt_spin == 1)then  
        do ! Run to equililibrium: cycles through the first 12 months or 365 days
          k = k + 1
          j = j + 1 
       
          if(k.eq.year_end+1)k = 1   ! 13 if monthly or 366 if daily
          if (test < 1E-8) exit
          YEAR = t_year(k)
          TEMP = t_tmp(k)
          RAIN = t_rain(k)
          PEVAP = t_evap(k)
         
          PC = t_PC(k)
          Pl_DPM_f = t_Pl_DPM_f(k)
          Pl_RPM_f = t_Pl_RPM_f(k)
         
          OA_DPM_f = t_OA_DPM_f(k)
          OA_RPM_f = t_OA_RPM_f(k)
          OA_Bio_f = t_OA_Bio_f(k)
          OA_Hum_f = t_OA_Hum_f(k)
         
         
          Pl_inp = t_Pl_inp(k)
          OA_inp = t_OA_inp(k)
         
          modernC = t_mod(k) / 100.0             
         
          call RothC(opt_tstep, opt_spin, DPM,RPM,Bio,Hum,IOM, SOC, 
     &      total_CO2, DPM_Rage, RPM_Rage, Bio_Rage, Hum_Rage, 
     &      Total_Rage, modernC, clay, depth,TEMP,RAIN,PEVAP,PC,
     &      Pl_DPM_f, Pl_RPM_f, OA_DPM_f, OA_RPM_f, OA_Bio_f, OA_Hum_f,
     &      Pl_inp, OA_inp, SMD, RM_TMP, RM_Moist, RM_PC, 
     &      opt_RMmoist, opt_SMDbare, silt, BD, OC, minRM_Moist)  
        
          if(mod(k, year_end)== 0)then 
            TOC0 = TOC1
            TOC1 = DPM+RPM+Bio+Hum
            test = abs(TOC1-TOC0)            
          endif    
         
        enddo
      
      else  ! if opt_spin is 2 or 3, set the dpm, rpm, bio, hum 
       dpm= dpm_init
       rpm= rpm_init
       bio= bio_init
       hum= hum_init
       soc= dpm + rpm + bio + hum + iom
       j = 1
          
      endif
       
      
      total_CO2 = 0.0 ! reset CO2 to zero after the equilibrium run
      
      write(91,9102) Year, j-1, DPM, RPM, Bio, Hum, IOM, SOC, total_CO2
     
9102     format(1x, i7, ',', i6, ',', 13x, ',', 15x, ',', 7x, ',',
     &         7x, ',', 8x, ',', 9x, ',', 7x, ',', 9x, ',',
     &           3x, ',', 6x, ',',f11.4, ',',f11.4,',', f11.4, ',',
     &        f11.4, ',',f11.4, ',',f11.4, ',',f11.4)   
     

C      
C run RothC to equilibrium: END
C
         
      Total_Delta = (exp(-Total_Rage/8035.0) - 1.0) * 1000.0   
      
      if(opt_spin==1)then
        write(71,102) year, j-1, DPM, RPM, Bio, Hum, iom, SOC,  
     &             total_CO2, Total_Delta
      else
        write(71,112) year, j-1, DPM, RPM, Bio, Hum, iom, SOC,  
     &             total_CO2
      endif
      
102   format(1x, i8, ',', 1x, i6, ',', 7(f11.4,','),  f8.2)     
112   format(1x, i8, ',', 1x, i6, ',', 6(f11.4,','),  f11.4)
C      
C run RothC for remaining timesteps to the end: START
C     
 !add an ifelse opt_spin = 2 here
 ! add an else here for opt_spin = 3
                                                             
      k_tstep = 0
      
      if(opt_spin==3)then
          start_loop = 1
      else
          start_loop= year_end+1
      endif
      
      do i = start_loop, nsteps, 1   ! 13 if monthly or 366 if daily
      
	  k_tstep = k_tstep + 1
            
        if(k_tstep.eq.year_end+1)k_tstep = 1   ! 13 if monthly or 366 if daily
         
        
         YEAR = t_year(i)
         TEMP = t_tmp(i)
         RAIN = t_rain(i)
         PEVAP = t_evap(i)
         
         PC = t_PC(i)
         Pl_DPM_f = t_Pl_DPM_f(i)
         Pl_RPM_f = t_Pl_RPM_f(i)
         
         OA_DPM_f = t_OA_DPM_f(i)
         OA_RPM_f = t_OA_RPM_f(i)
         OA_Bio_f = t_OA_Bio_f(i)
         OA_Hum_f = t_OA_Hum_f(i)
         
         Pl_inp = t_Pl_inp(i)
         OA_Inp = t_OA_inp(i)
         
         modernC = t_mod(i) / 100.0
           
         call RothC(opt_tstep, opt_spin, DPM,RPM,Bio,Hum,IOM, SOC, 
     &     total_CO2, DPM_Rage, RPM_Rage, Bio_Rage, Hum_Rage, 
     &     Total_Rage, modernC, clay, depth,TEMP,RAIN,PEVAP,PC,
     &     Pl_DPM_f, Pl_RPM_f, OA_DPM_f, OA_RPM_f, OA_Bio_f, OA_Hum_f,
     &     Pl_inp, OA_inp, SMD, RM_TMP, RM_Moist, RM_PC, 
     &     opt_RMmoist, opt_SMDbare, silt, BD, OC, minRM_Moist)    
         
         Total_Delta = (exp(-Total_Rage/8035.0) - 1.0) * 1000.0
         
         
         write(91,9103) Year, k_tstep, Pl_inp, OA_inp, TEMP,RM_TMP, 
     &        RAIN, PEVAP, SMD, RM_Moist, PC, RM_PC,
     &        DPM,RPM,Bio,Hum, IOM, SOC, total_CO2
     
9103     format(1x, i7, ',', i6, ',', f13.3, ',',f15.3, ',',f7.1, ',',
     &         f7.4, ',',f8.1, ',',f9.1, ',',f7.2, ',', f9.4, ',',
     &           i3, ',',f6.1, ',',f11.4, ',',f11.4,',', f11.4, ',',
     &        f11.4, ',',f11.4, ',',f11.4, ',',f11.4)        


      if(mod(i, year_end)== 0)then     ! print out results once a year
        if(opt_spin==1)then
          write(71,103) year, DPM, RPM, Bio, Hum, IOM, SOC, total_CO2, 
     &                Total_Delta
        else
          write(71,113) year, DPM, RPM, Bio, Hum, IOM, SOC, total_CO2
        endif
      endif
         
      enddo  
C      
C run RothC for remaining timesteps to the end: END
C   
      
      call cpu_time (time_end)
      
!      write(81,*) 'Time of operation was ', 
!     $    time_end - time_begin, ' seconds'    
           
 103  format(1x,  i7, ',', 6x, '12,', 7(f11.4, ','), f8.2)  
 113  format(1x,  i7, ',', 6x, '12,', 6(f11.4, ','), f11.4)  
 
      close (71)
      close (91)
  
      stop
      
      end
      
            