C******************************************************************************
C  RothC model
C
C  February 2024
C
C  
C  Kevin Coleman
C
C  This is the code for RothC 
C
C  June 2025 this is a branch for the Farina (2013) version of the model
C
C INPUTS: 
C
C clay:        clay content of the soil (units: %)
C depth:       depth of topsoil (units: cm)
C IOM:         inert organic matter (t C /ha)
C nsteps:      number of timesteps 
C
C The following are needed for the Farina modification to the model (Farina et al, 2013, Geoderma. 200, 18-30, 10.1016/j.geoderma.2013.01.021)
C
C slit:        silt content of the soil (units: %) 
C BD:          bulk density (units: g/cm3) for the farina (2013) version
C OC:          organic carbon (units: %) for the farina (2013) version
C minRM_Moist: the minimum value the rate modifying factor for moisture can be (units: -, default=0.2)
C
C The following switches are needed to allow the user to choose which model option to run
C
C opt_RMmoist !  1: Standard RothC soil water parameters,
C             !  2: Van Genuchten soil properties and soil is allowed to be drier (ie hygroscopic / capillary water, -1000bar)
C             !  3: Van Genuchten soil properties, but uses the Standard RothC soil water function
C      
C opt_SMDbare !  1: Standard RothC bareSMD, 
C             !  2: bareSMD is set to wilting point -15bar (could be better for dry soils)
C
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
C  SMD:       soil moisture deficit (mm per soil depth)  SMD
C  RM_TMP:    rate modifying fator for temperature (0.0 - ~5.0)
C  RM_Moist:  rate modifying fator for moisture (0.0 - 1.0)
C  RM_PC:     rate modifying fator for plant retainment (0.6 or 1.0)
C
C******************************************************************************      
      Subroutine RothC(timeFact, DPM,RPM,BIO,HUM,IOM, SOC, total_CO2, 
     &         DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage,
     &         modernC, clay, depth,TEMP,RAIN,PEVAP,PC,DPM_RPM,
     &         C_Inp, FYM_Inp, SMD, RM_TMP, RM_Moist, RM_PC, 
     &         opt_RMmoist, opt_SMDbare, silt, BD, OC, minRM_Moist) 
         
      
      
      implicit none
      
      integer timeFact
      
      integer PC
      
      integer opt_RMmoist !  1: Standard RothC soil water parameters,
                          !  2: Van Genuchten soil properties and soil is allowed to be drier (ie hygroscopic / capillary water, -1000bar)
                          !  3: Van Genuchten soil properties, but uses the Standard RothC soil water function
      
      integer opt_SMDbare !  1: Standard RothC bareSMD, 
                          !  2: bareSMD is set to wilting point -15bar (could be better for dry soils)

      real*8 DPM,RPM,BIO,HUM,IOM,SOC, total_CO2
      
      real*8 DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage
      
      real*8 modernC
      
      real*8 clay, depth
      
      real*8 silt  ! silt content (units: %) needed for the farina (2013) version
     
      real*8 BD  ! bulk density (units: g/cm3) needed for the farina (2013) version
               
      real*8 OC  ! organic carbon (units: %) needed for the farina (2013) version
      
      real*8 minRM_Moist  ! (units: -, default=0.2) needed for the farina (2013) version
      
      real*8 DPM_RPM
      
      real*8 TEMP, RAIN, PEVAP
      
      real*8 C_Inp, FYM_Inp
      
      real*8 SMD
       
      real*8 RM_TMP, RM_Moist, RM_PC   
      
      real*8 RateM
      
       
C Calc RMF's i've create a subroutine for each to make it easier to add others or replace existing with new.     
      call RMF_Tmp(TEMP, RM_TMP)
      call RMF_Moist(RAIN, PEVAP, clay, depth, PC, SMD, RM_Moist, 
     &              opt_RMmoist, opt_SMDbare, silt, BD, OC, minRM_Moist)
      call RMF_PC(PC,RM_PC)

      
      
C combine RMF's into one.      
      RateM = RM_TMP*RM_Moist*RM_PC
      

      
      call decomp(timeFact, DPM,RPM,BIO,HUM, IOM, SOC, total_CO2,   
     &     DPM_Rage, RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage, modernC,
     &     RateM, clay, C_Inp, FYM_Inp, DPM_RPM)
     
     
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
      Subroutine RMF_Moist (RAIN, PEVAP, clay, depth, PC, SMD, RM_Moist,
     &         opt_RMmoist, opt_SMDbare, silt, BD, OC, minRM_Moist)
      
      implicit none
      
      integer opt_RMmoist !  1: Standard RothC soil water parameters,
                          !  2: Van Genuchten soil properties and soil is allowed to be drier (ie hygroscopic / capillary water, -1000bar)
                          !  3: Van Genuchten soil properties, but uses the Standard RothC soil water function
      
      integer opt_SMDbare !  1: Standard RothC bareSMD, 
                          !  2: bareSMD is set to wilting point -15bar (could be better for dry soils)
      
      integer PC

      real*8 RAIN, PEVAP
      real*8 clay, depth
      
      real*8 silt, BD, OC  ! needed for the farina (2013) version
      
      real*8 RM_Moist
      
      real*8 SMD
      
      real*8 DF
      
      real*8 SMDMax, SMDMaxAdj, SMD1bar, SMDBare 
      
      real*8 SMD15bar, SMD15barAdj, SMD1000bar  ! for the farina (2013) version
      
      real*8, parameter :: RMFMax = 1.0
      
      real*8 minRM_Moist

      real*8 X0, X1, X2, X3, D2

       
!C calc soil water functions properties       
      IF(opt_RMmoist.eq.1)THEN                   ! 1: use standard RothC soil water properties 
        SMD15bar=-(20+1.3*clay-0.01*(clay*clay))
        SMD15barAdj = SMD15bar * depth / 23.0
        SMD1bar = 0.444 * SMD15barAdj
      ELSE									   ! 2 or 3: Van Genuchten soil properties  !!  and soil is allowed to be drier 
        CALL CALC_SM_VG(clay, silt, BD, OC, depth, 
     &                       X0, X1, X2, X3)
        SMD15bar = X2  ! X2 has been adjusted for depth in CALC_SM_SM
        SMD15barAdj = SMD15bar
        SMD1bar = X1
        SMD1000bar = X3
      ENDIF 
      
      IF(opt_SMDbare.EQ.1)THEN                   ! Standard RothC bareSMD
         IF(opt_RMmoist.EQ.1)THEN       ! 
           SMDBare = 0.556 * SMD15barAdj
	   ELSE  ! i think we need this option if option 3 standard  
           SMDBare = SMD15barAdj - (0.6388/0.8) * (SMD15barAdj-SMD1bar)
	   ENDIF
	ELSE                                       ! bareSMD is set to wilting point -15bar (could be better for dry soils)
        SMDBare = SMD15barAdj
	ENDIF    
      
      IF(opt_RMmoist.EQ.1.or.opt_RMmoist.EQ.3)THEN                   ! use standard RothC soil water properties   
        SMDMaxAdj = SMD15barAdj
      ELSEIF(opt_RMmoist.EQ.2)THEN               ! Van Genuchten soil properties and soil is allowed to be drier 
        SMDMaxAdj = SMD1000bar
      ENDIF      
      
      DF = RAIN - 0.75 * PEVAP
      
      IF(PC.eq.1)THEN
        SMD = MAX(SMDMAXadj,MIN(0.0,SMD+DF))
      ELSE
        SMD = MAX(MIN(SMDBare,SMD),MIN(0.0,SMD+DF))
      ENDIF
        
      IF(opt_RMmoist.eq.1.or.opt_RMmoist.eq.3)THEN                   ! use standard RothC soil water properties 
      
          IF(SMD.gt.SMD1bar)THEN
            RM_Moist = 1.0
          ELSE
            RM_Moist = (minRM_Moist + (RMFMax - minRM_Moist) * 
     &                (SMD15barAdj - SMD) / (SMD15barAdj - SMD1bar) ) 
          ENDIF
          
      ELSE                                      ! Van Genuchten soil properties and soil is allowed to be drier 
          IF(SMD.gt.SMD1bar)THEN
            RM_Moist = 1.0
          ELSEIF(SMD.gt.SMD15barAdj)THEN
            RM_Moist = (minRM_Moist + (RMFMax - minRM_Moist) * 
     &                (SMD15barAdj - SMD) / (SMD15barAdj - SMD1bar) )   
          ELSE
            RM_Moist = minRM_Moist       
          ENDIF
      ENDIF

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
      Subroutine decomp(timeFact, DPM,RPM,BIO,HUM, IOM, SOC, total_CO2,
     &       DPM_Rage,RPM_Rage, Bio_Rage, HUM_Rage, Total_Rage, modernC,
     &       RateM, clay, C_Inp, FYM_Inp, DPM_RPM)
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
      
      real*8 DPM,RPM,BIO,HUM,IOM, SOC, total_CO2
      
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
      
      total_CO2 = total_CO2 + DPM_co2 + RPM_co2 + Bio_co2 + Hum_co2
      
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
      
C
C ******************************************************************************
C ******************************************************************************
C
C
      SUBROUTINE CALC_SM_VG(clayper, siltper, BD, OC, depth, 
     &                       X0, X1, X2, X3)

      IMPLICIT none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculates the soil properties using Wosten et al. (1999) based on van Genuchten (1980).
c
c van Genuchten, M.T., 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. 
c SSSAJ, 44, 892-898, 
c
c Wosten, J.H.M., Lilly, A., Nemes, A., Le Bas, C., 1999. 
c Development and use of a database of hydraulic properties of European soils. 
c Geoderma. 90, 169-185, http://dx.doi.org/10.1016/s0016-7061(98)00132-3
c
c
c    clay - clay (%)
c    silt - silt (%)
c    BD  - Bulk density (g cm-3)
c    OC - Organic carbon (%)
c    depth - depth of soil sample (cm) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision clayper, siltper, BD, OC, depth

	integer mbars(5)

	integer i 

	double precision alpha, thetaS, n, thetaR, m, ksat, l_star, l, t

      double precision wc(5)

      double precision wcSAT
	double precision wcFC
	double precision wc1
	double precision wcWP
	double precision wc1000

c     X0 - Soil moisture deficit (between FC and Saturation) this is positive, so that FC reminds zero
c     X1 - Soil moisture deficit (between FC and 1 bar)
c     X2 - maximum Soil moisture deficit (between FC and WP)
c     X3 - Soil moisture deficit (between FC and 1000 bar)
      double precision X0, X1, X2, X3

      mbars(1) = 0 
      mbars(2) = 50
	mbars(3) = 1000
	mbars(4) = 15000
	mbars(5) = 1000000

	t = 1   ! note t = topsoil, in Wosten et al 1999, topsoil and subsoil have values 1 or 0
              ! RothC only models topsoils so t is always 1

      alpha=EXP(-14.96+0.03135*clayper+0.0351*siltper+0.646*(OC*1.72)
     &     +15.29*BD-0.192*t-4.671*BD**2-0.000781*clayper**2
     &     -0.00687*(OC*1.72)**2
     &     +0.0449*(OC*1.72)**-1+0.0663*log(siltper)
     &     +0.1482*log(OC*1.72)
     &     -0.04546*BD*siltper-0.4852*BD*(OC*1.72)+0.00673*clayper*t)

      thetaS=(0.7919+0.001691*clayper-0.29619*BD-0.000001491*siltper**2
     &     +0.0000821*(OC*1.72)**2+0.02427*clayper**-1
     &     +0.01113*siltper**-1+0.01472*log(siltper)
     &     -0.0000733*(OC*1.72)*clayper-0.000619*BD*clayper
     &     -0.001183*BD*(OC*1.72)-0.0001664*siltper*t)

      n=EXP(-25.23 -0.02195*clayper +0.0074*siltper -0.194*(OC*1.72)
     &     +45.5*BD-7.24*BD**2 +0.0003658*clayper**2
     &     +0.002885*(OC*1.72)**2 -12.81*BD**-1 -0.1524*siltper**-1
     &     -0.01958*(OC*1.72)**-1 -0.2876*log(siltper)
     &     -0.0709*log(OC*1.72) -44.6*log(BD) -0.02264*BD*clayper
     &     +0.0896*BD*(OC*1.72) +0.00718*clayper*t)+1

      thetaR=0.01

      m=1-1/n

      ksat=EXP(7.755 +0.0352*siltper +0.93*t -0.967*BD**2 
     &     -0.000484*clayper**2 -0.000322*siltper**2
     &     +0.001*siltper**-1 -0.0748*(OC*1.72)**-1
     &     -0.643*log(siltper) -0.01398*BD*clayper -0.1673*BD*(OC*1.72)
     &     +0.02986*clayper*t -0.03305*siltper*t)

      l_star=(0.0202 +0.0006193*clayper**2 -0.001136*(OC*1.72)**2
     &     -0.2316*log(OC*1.72) -0.03544*BD*clayper +0.00283*BD*siltper
     &     +0.0488*BD*OC*1.72)

      l=10*(EXP(l_star)-1)/(EXP(l_star)+1)


      do i = 1, 5
        wc(i)=thetaR + (thetaS-thetaR)/ (1+(alpha*mbars(i))**n)**m
      enddo

      wcSAT  = wc(1)
	wcFC   = wc(2)
	wc1    = wc(3)
	wcWP   = wc(4)
	wc1000 = wc(5)
     
      X0 = (wcSAT  - wcFC)*10*depth
	X1 = (wc1    - wcFC)*10*depth
	X2 = (wcWP   - wcFC)*10*depth
	X3 = (wc1000 - wcFC)*10*depth
      
	RETURN
      END
