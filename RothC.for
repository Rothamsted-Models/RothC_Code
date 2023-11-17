C******************************************************************************
C  RothC model
C
C  December 2021

C  
C  Kevin Coleman
C
C  This is the code for RothC 
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
C
C OUTPUTS:

C All pools are carbon and not organic matter 
C
C  DPM:   Decomposable Plant Material (units: t C /ha)
C  RPM:   Resistant Plant Material    (units: t C /ha)
C  Bio:   Microbial Biomass           (units: t C /ha)
C  Hum:   Humified Organic Matter     (units: t C /ha)
C  IOM:   Inert Organic Matter        (units: t C /ha)
C  SOC:   Soil Organic Matter / Total organic Matter (units: t C / ha)

C  DPM_Rage:   radiocarbon age of DPM
C  RPM_Rage:   radiocarbon age of RPM
C  Bio_Rage:   radiocarbon age of Bio
C  HUM_Rage:   radiocarbon age of Hum
C  Total_Rage: radiocarbon age of SOC (/ TOC)
C
C  SWC:       soil moisture deficit (mm per soil depth)
C  RM_TMP:    rate modifying fator for temperature (0.0 - ~5.0)
C  RM_Moist:  rate modifying fator for moisture (0.0 - 1.0)
C  RM_PC:     rate modifying fator for plant retainment (0.6 or 1.0)

C
C******************************************************************************      
      Subroutine RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC,
     &         DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage,
     &         modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM,
     &         C_Inp, FYM_Inp, SWC, RM_TMP, RM_Moist, RM_PC) 
         
      
      
      implicit none
      
      integer timeFact
      
      integer PC
      
      real*8 DPM,RPM,BIO,HUM,IOM,SOC
      
      real*8 DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage
      
      real*8 modernC
      
      real*8 clay, depth
      
      real*8 DPM_RPM
      
      real*8 TEMP, RAIN, PEVAP
      
      real*8 C_Inp, FYM_Inp
      
      real*8 SWC
       
      real*8 RM_TMP, RM_Moist, RM_PC   
      
      real*8 RateM
      
       
C Calc RMF's i've create a subroutine for each to make it easier to add others or replace existing with new.     
      call RMF_Tmp(TEMP, RM_TMP)
      call RMF_Moist(RAIN, PEVAP, clay, depth, PC, SWC, RM_Moist)
      call RMF_PC(PC,RM_PC)

C combine RMF's into one.      
      RateM = RM_TMP*RM_Moist*RM_PC
      

      
      call decomp(timeFact, DPM,RPM,BIO,HUM, IOM, SOC, DPM_Rage,  
     &       RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage, modernC, RateM,
     &       clay, C_Inp, FYM_Inp, DPM_RPM)
     
     
      return
      end
C
C**********************************************************************
C    Calculates the temperature modifying factor   
C**********************************************************************     
C
      Subroutine RMF_Tmp (TEMP, RM_TMP)
C
      implicit none

      real*8 TEMP
      real*8 RM_TMP
    
      
      IF(TEMP.LT.-5.0)THEN
        RM_TMP=0.0
      ELSE
        RM_TMP=47.91/(EXP(106.06/(TEMP+18.27))+1.0)
      END IF

      RETURN
      END
C
C**********************************************************************
C       Calculates the moisture modifying factor  
C**********************************************************************     
C
      Subroutine RMF_Moist (RAIN, PEVAP, clay, depth, PC, SWC, RM_Moist)

      implicit none
      
      integer PC

      real*8 RAIN, PEVAP
      real*8 clay, depth
      real*8 RM_Moist
      
      real*8 SWC
      
      real*8 DF
      
      real*8 SMDMax, SMDMaxAdj, SMD1bar, SMDBare 
      
      real*8, parameter :: RMFMax = 1.0, RMFMin = 0.2

C calc soil water functions properties
      SMDMax=-(20+1.3*clay-0.01*(clay*clay))
      SMDMaxAdj = SMDMax * depth / 23.0
      SMD1bar = 0.444 * SMDMaxAdj
      SMDBare = 0.556 * SMDMaxAdj
      
      DF = RAIN - 0.75 * PEVAP
      
      if(PC.eq.1)then
        SWC = MAX(SMDMAXadj,MIN(0.0,SWC+DF))
      else
        SWC = MAX(MIN(SMDBare,SWC),MIN(0.0,SWC+DF))
      endif  
      
      if(SWC.gt.SMD1bar)then
        RM_Moist = 1.0
      else
        RM_Moist = (RMFMin + (RMFMax - RMFMin) * 
     &              (SMDMaxAdj - SWC) / (SMDMaxAdj - SMD1bar) )   
      endif

      RETURN
      END
C
C**********************************************************************
C      Calculates the plant retainment modifying factor  
C**********************************************************************     
C      
      Subroutine RMF_PC (PC, RM_PC)
C
      implicit none

      integer PC
      real*8 RM_PC
      
      IF(PC.eq.0)THEN
        RM_PC = 1.0
      ELSE
        RM_PC = 0.6
      END IF

      RETURN
      END
C
C**********************************************************************
C      calculates the decomposition and radiocarbon 
C**********************************************************************     
C
      Subroutine decomp(timeFact, DPM,RPM,BIO,HUM, IOM, SOC, DPM_Rage, 
     &       RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage, modernC, RateM,  
     &       clay, C_Inp, FYM_Inp, DPM_RPM)
C
      implicit none
      
      real*8, parameter :: zero = 0e-8
C rate constant are params so don't need to be passed
!      real*8, parameter :: tstep = 1/12.0
      real*8, parameter :: DPM_k = 10.0,  RPM_k = 0.3
      real*8, parameter :: Bio_k = 0.66,  Hum_k = 0.02 
    
      real*8, parameter ::  conr = log(2.0) / 5568.0
!     real*8, parameter ::  exc = exp(-conr*tstep) 
      
      integer timeFact
      
      real*8 tstep, exc
      
      real*8 DPM,RPM,BIO,HUM,IOM, SOC
      
      real*8 modernC
      
      real*8 RateM
      real*8 Clay
      
      real*8 C_Inp, FYM_Inp, DPM_RPM
      
      real*8 PI_C_DPM, PI_C_RPM   
      real*8 FYM_C_DPM, FYM_C_RPM, FYM_C_Hum

C C that remains in each pool       
      real*8 DPM1,RPM1,BIO1,HUM1 
           
C C amount that will become CO2, Bio and Hum (_d = delta change)
      real*8 DPM_d,RPM_d,BIO_d,HUM_d      
      
c _co2: amount of pool that becomes CO2 
      real*8 DPM_co2, RPM_co2, Bio_co2, Hum_co2
      
c _bio: amount of pool that becomes bio 
      real*8 DPM_bio, RPM_bio, Bio_bio, Hum_bio
      
c _hum amount of pool that becomes hum 
      real*8 DPM_hum, RPM_hum, Bio_hum, Hum_hum 
      
      real*8 DPM_Rage,RPM_Rage,Bio_Rage,HUM_Rage,IOM_Rage,Total_Rage
      
      real*8 DPM_Ract,RPM_Ract,Bio_Ract,Hum_Ract,IOM_Ract,Total_Ract
      
      real*8 PI_DPM_Ract, PI_RPM_Ract
      
      real*8 FYM_DPM_Ract, FYM_RPM_Ract, FYM_Hum_Ract      
      
      real*8 DPM_Ract_new, RPM_Ract_new, Bio_Ract_new, Hum_Ract_new
      
      real*8 DPM_Bio_Ract, RPM_Bio_Ract, Bio_Bio_Ract, Hum_Bio_Ract 
      real*8 DPM_Hum_Ract, RPM_Hum_Ract, Bio_Hum_Ract, Hum_Hum_Ract
 
      real*8 X
 
      tstep = 1.0/timeFact    ! monthly 1/12, or daily 1/365  
      
      exc = exp(-conr*tstep) 
      
      IOM_Rage =50000.0
      
 
C C decomposition
      DPM1 = DPM * exp(-RateM*DPM_k*tstep)
      RPM1 = RPM * exp(-RateM*RPM_k*tstep)      
      Bio1 = Bio * exp(-RateM*Bio_k*tstep)      
      Hum1 = Hum * exp(-RateM*Hum_k*tstep) 

      
      DPM_d = DPM - DPM1
      RPM_d = RPM - RPM1      
      Bio_d = Bio - Bio1
      Hum_d = Hum - Hum1 
      

      X=1.67*(1.85+1.60*EXP(-0.0786*Clay))

 
C proportion C from each pool into CO2, Bio and Hum      
      DPM_co2 = DPM_d * (x / (x+1))
      DPM_bio = DPM_d * (0.46 / (x+1))
      DPM_hum = DPM_d * (0.54 / (x+1))
      
      RPM_co2 = RPM_d * (x / (x+1))
      RPM_bio = RPM_d * (0.46 / (x+1))
      RPM_hum = RPM_d * (0.54 / (x+1))    
      
      Bio_co2 = Bio_d * (x / (x+1))
      Bio_bio = Bio_d* (0.46 / (x+1))
      Bio_hum = Bio_d * (0.54 / (x+1))
      
      Hum_co2 = Hum_d * (x / (x+1))
      Hum_bio = Hum_d * (0.46 / (x+1))
      Hum_hum = Hum_d * (0.54 / (x+1))  
           
      
C update C pools  
      DPM = DPM1
      RPM = RPM1
      Bio = Bio1 + DPM_bio + RPM_bio + Bio_bio + Hum_bio
      Hum = Hum1 + DPM_hum + RPM_hum + Bio_hum + Hum_hum    
      
C split plant C to DPM and RPM 
      PI_C_DPM = DPM_RPM / (DPM_RPM + 1.0) * C_Inp
      PI_C_RPM =     1.0 / (DPM_RPM + 1.0) * C_Inp

C split FYM C to DPM, RPM and Hum 
      FYM_C_DPM = 0.49*FYM_Inp
      FYM_C_RPM = 0.49*FYM_Inp      
      FYM_C_Hum = 0.02*FYM_Inp   
      
C add Plant C and FYM_C to DPM, RPM and Hum   
      DPM = DPM + PI_C_DPM + FYM_C_DPM
      RPM = RPM + PI_C_RPM + FYM_C_RPM  
      Hum = Hum + FYM_C_Hum
      
C calc new ract of each pool      
      DPM_Ract = DPM1 *exp(-conr*DPM_Rage)
      RPM_Ract = RPM1 *exp(-conr*RPM_Rage) 
      
      Bio_Ract = Bio1 *exp(-conr*Bio_Rage)
      DPM_Bio_Ract = DPM_Bio * exp(-conr*DPM_Rage)
      RPM_Bio_Ract = RPM_Bio * exp(-conr*RPM_Rage)
      Bio_Bio_Ract = Bio_Bio * exp(-conr*Bio_Rage)
      Hum_Bio_Ract = Hum_Bio * exp(-conr*Hum_Rage)
      
      Hum_Ract = Hum1 *exp(-conr*Hum_Rage)   
      DPM_Hum_Ract = DPM_Hum * exp(-conr*DPM_Rage)
      RPM_Hum_Ract = RPM_Hum * exp(-conr*RPM_Rage)
      Bio_Hum_Ract = Bio_Hum * exp(-conr*Bio_Rage)
      Hum_Hum_Ract = Hum_Hum * exp(-conr*Hum_Rage)
      
      IOM_Ract = IOM *exp(-conr*IOM_Rage) 
      
C assign new C from plant and FYM the correct age   
      PI_DPM_Ract = modernC * PI_C_DPM
      PI_RPM_Ract = modernC * PI_C_RPM
      
      FYM_DPM_Ract = modernC * FYM_C_DPM
      FYM_RPM_Ract = modernC * FYM_C_RPM
      FYM_Hum_Ract = modernC * FYM_C_Hum          
      
C update ract for each pool        
      DPM_Ract_new = FYM_DPM_Ract + PI_DPM_Ract + DPM_Ract*exc
      RPM_Ract_new = FYM_RPM_Ract + PI_RPM_Ract + RPM_Ract*exc    
      
      Bio_Ract_new = (Bio_Ract + DPM_Bio_Ract + RPM_Bio_Ract + 
     &                Bio_Bio_Ract + Hum_Bio_Ract )*exc
      
      Hum_Ract_new = FYM_Hum_Ract + (Hum_Ract + DPM_Hum_Ract +
     &               RPM_Hum_Ract + Bio_Hum_Ract + Hum_Hum_Ract )*exc  
      
      
      SOC = DPM + RPM + Bio + Hum + IOM      
      
      Total_Ract = DPM_RACT_new + RPM_Ract_new +
     &           Bio_Ract_new + HUM_Ract_new + IOM_Ract
      

C calculate rage of each pool.      
      if(DPM.le.zero)then
        DPM_Rage = zero
      else
        DPM_Rage = ( log(DPM/DPM_Ract_new) ) / conr
      endif
      
      if(RPM.le.zero)then
        RPM_Rage = zero
      else
        RPM_Rage = ( log(RPM/RPM_Ract_new) ) / conr 
      endif
      
      if(Bio.le.zero)then
        Bio_Rage = zero
      else
        Bio_Rage = ( log(Bio/Bio_Ract_new) ) / conr
      endif
      
      if(Hum.le.zero)then
        Hum_Rage = zero
      else
        Hum_Rage = ( log(Hum/Hum_Ract_new) ) / conr
      endif
        
      if(SOC.le.zero)then
        Total_Rage = zero
      else
        Total_Rage = ( log(SOC/Total_Ract) ) / conr    
      endif
      
      RETURN
      END