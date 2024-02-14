C******************************************************************************
C  Wrapper for RothC model 
C
C  November 2023
C
C 
C  
C  Kevin Coleman
C
C******************************************************************************   
C
C INPUTS: 
C
C clay:  clay content of the soil (units: %)
C depth: depth of topsoil (units: cm)
C IOM: inert organic matter (t C /ha)
C nsteps: number of timesteps 
C
C year:    year
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
      
      real*8 DPM, RPM, BIO, HUM, IOM, SOC
      
      real*8 DPM_Rage,RPM_Rage,BIO_Rage,HUM_Rage,IOM_Rage,Total_Rage
      
      
      real*8 DPM_Delta, RPM_Delta, Bio_Delta, Hum_Delta, IOM_Delta
      real*8 Total_Delta
      
      integer YEAR
      
      real*8 TEMP, RAIN, PEVAP
      
      real*8 RM_TMP, RM_Moist, RM_PC
      
      real*8 modernC
      
      ! C_inp (Carbon input to the soil each month units: t C /ha)
      
      real*8 C_inp, FYM_Inp, DPM_RPM  
      
      real*8 SWC      
      
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
      
      DPM_Rage = 0.0
      RPM_Rage = 0.0
      BIO_Rage = 0.0
      Hum_Rage = 0.0
      IOM_Rage = 50000.0  

C set initial soil water content (deficit) 
      SWC = 0.0
C
C READ IN INPUT DATA: START
C
C read in RothC input data file: data will be passed from other programs at some point  
      open(11, file='RothC_input.dat', status='unknown')    
!      open(11, file='RothC_input_Morocco.dat', status='unknown') 
	read(11,*)
	read(11,*)      
	read(11,*) 
	read(11,*)             
      read(11,*)clay, depth, iom, nsteps
      read(11,*)
      read(11,*)  

	do i = 1, nsteps
	  read(11,*)t_year(i), t_mod(i), t_tmp(i),t_rain(i),t_evap(i),
     &            t_C_Inp(i), t_FYM_Inp(i), t_PC(i), t_DPM_RPM(i)
      enddo
      
      close(11)   
      
      
      open(71, file='year_results.dat', status ='unknown')     
      open(91, file='month_results.dat', status ='unknown')
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
7100  format(5x, 'Year', 3x,  'Month ', 
     &  1x, 'DPM_t_C_ha ', 1x,  'RPM_t_C_ha ', 
     &  1x, 'BIO_t_C_ha ', 1x, 'HUM_t_C_ha ',
     &  1x, 'IOM_t_C_ha ', 1x, 'SOC_t_C_ha ', 
     &  1x, 'deltaC')     
             
      write(71,101) j, DPM, RPM, Bio, Hum, iom, SOC
101   format(1x, '       0', 1x, i7, 6f12.4, ' -998.02')  

      timeFact = 12 ! monthly
          
      test = 100.0
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
         
         call RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC, 
     &      DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage, 
     &      modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM,
     &      C_Inp, FYM_Inp, SWC, RM_TMP, RM_Moist, RM_PC)  
        
         if(mod(k, timeFact)== 0)then 
           TOC0 = TOC1
           TOC1 =DPM+RPM+Bio+Hum
           test = abs(TOC1-TOC0)            
         endif       
         
      enddo
C      
C run RothC to equilibrium: END
C
         
      Total_Delta = (exp(-Total_Rage/8035.0) - 1.0) * 1000.0   
      
      write(71,102) year, j-1, DPM, RPM, Bio, Hum, iom, SOC, Total_Delta
      
102   format(1x, i8, 1x, i7, 6f12.4 , f8.2)             

       write(91,9100)
9100  format(4x, 'Year ',1x,  'Month ',1x, 'C_Inp_t_C_ha ', 
     &  1x,  'FYM_Inp_t_C_ha ', 1x,  'TEMP_C ', 1x, 'RM_TMP ',
     &  1x, 'RAIN_mm ', 1x, 'PEVAP_mm ',1x, 'SWC_mm ',
     &  1x,'RM_Moist ', 1x, 'PC ', 1x,  'RM_PC ',  
     &  1x, 'DPM_t_C_ha ', 1x,  'RPM_t_C_ha ', 
     &  1x, 'BIO_t_C_ha ', 1x, 'HUM_t_C_ha ',
     &  1x, 'IOM_t_C_ha ', 1x, 'SOC_t_C_ha ')
                                                             
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
           
         call RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC, 
     &      DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage, 
     &      modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM,
     &      C_Inp, FYM_Inp, SWC, RM_TMP, RM_Moist, RM_PC )  
         
         Total_Delta = (exp(-Total_Rage/8035.0) - 1.0) * 1000.0
         
         
         write(91,9101) Year, k_month, C_Inp, FYM_Inp, TEMP,RM_TMP, 
     &        RAIN, PEVAP, SWC, RM_Moist, PC, RM_PC,
     &        DPM,RPM,BIO,HUM, IOM, SOC
     
9101     format(1x, 2i7, f14.3, f16.3, f8.1, f8.4, f9.1, f10.1, f8.2, 
     &              f10.4, i4, f7.1, f12.4, f12.4, f12.4, f12.4, f12.4,
     &              f12.4)         


      if(mod(i, timeFact)== 0)then     ! print out results once a year
        write(71,103) year, DPM, RPM, Bio, Hum, IOM, SOC, Total_Delta
      endif
         
      enddo  
C      
C run RothC for months 13 to the end: END
C   
      
      call cpu_time (time_end)
      
      write(81,*) 'Time of operation was ', 
     $    time_end - time_begin, ' seconds'
           
 103  format(1x,  i8, 6x, '12', 6f12.4, f8.2)  
 
      close (71)
      close (91)
  
      stop
      
      end
      
            