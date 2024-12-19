      ! POSTPROCESS W2 RESULTS

      PROGRAM TEMPERATURE

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      !===================================================================================================================
      ! 223200 = 24*31*12*25, 24 HRS, MAX 31 DAYS, 12 MONTHS, 25 YEARS (WATER YEAR 1979~2003)
	! 775 = 31*25

c      real*8 ttest,ttest1

      REAL*8 TEMP_ML(4,223200)  ! temperature
     +	   ,Q_ML(4,223200)       ! flow rate
     +       ,V_ML(3,223200)       ! reservoir volume
     +       ,E_ML(4,223200)       ! elevation of structure

      ! _1: DAILY AVERAGE RESULTS, _2: MAXIMUM HOURLY RESULTS FOR EACH DAY,  _3: MINIMUM HOURLY RESULTS FOR EACH DAY
 
	REAL*8 TOCT_ML_1(4,775),TNOV_ML_1(4,775),TDEC_ML_1(4,775),
     +       TJAN_ML_1(4,775),TFEB_ML_1(4,775),TMAR_ML_1(4,775),
     +       TAPR_ML_1(4,775),TMAY_ML_1(4,775),TJUN_ML_1(4,775),
     +       TJUL_ML_1(4,775),TAUG_ML_1(4,775),TSEP_ML_1(4,775)   

	REAL*8 TOCT_ML_2(4,775),TNOV_ML_2(4,775),TDEC_ML_2(4,775),
     +       TJAN_ML_2(4,775),TFEB_ML_2(4,775),TMAR_ML_2(4,775),
     +       TAPR_ML_2(4,775),TMAY_ML_2(4,775),TJUN_ML_2(4,775),
     +       TJUL_ML_2(4,775),TAUG_ML_2(4,775),TSEP_ML_2(4,775)  

	REAL*8 TOCT_ML_3(4,775),TNOV_ML_3(4,775),TDEC_ML_3(4,775),
     +       TJAN_ML_3(4,775),TFEB_ML_3(4,775),TMAR_ML_3(4,775),
     +       TAPR_ML_3(4,775),TMAY_ML_3(4,775),TJUN_ML_3(4,775),
     +       TJUL_ML_3(4,775),TAUG_ML_3(4,775),TSEP_ML_3(4,775)    

      REAL*8 QOCT_ML_1(4,775),QNOV_ML_1(4,775),QDEC_ML_1(4,775),
     +       QJAN_ML_1(4,775),QFEB_ML_1(4,775),QMAR_ML_1(4,775),
     +       QAPR_ML_1(4,775),QMAY_ML_1(4,775),QJUN_ML_1(4,775),
     +       QJUL_ML_1(4,775),QAUG_ML_1(4,775),QSEP_ML_1(4,775)            

      REAL*8 EOCT_ML_1(4,775),ENOV_ML_1(4,775),EDEC_ML_1(4,775),
     +       EJAN_ML_1(4,775),EFEB_ML_1(4,775),EMAR_ML_1(4,775),
     +       EAPR_ML_1(4,775),EMAY_ML_1(4,775),EJUN_ML_1(4,775),
     +       EJUL_ML_1(4,775),EAUG_ML_1(4,775),ESEP_ML_1(4,775)

      OPEN (UNIT=2,FILE='..\..\str_br1.opt',STATUS='OLD')
      OPEN (UNIT=4,FILE='..\..\Volume_wb1.opt',STATUS='OLD')

      OPEN (UNIT=21,FILE='str_br1_1_ave.opt')	! output average hourly results for each day, ML spillway
      OPEN (UNIT=22,FILE='str_br1_2_ave.opt') ! output average hourly results for each day, F-K Canal
	OPEN (UNIT=23,FILE='str_br1_3_ave.opt') ! output average hourly results for each day, Madera Canal
	OPEN (UNIT=24,FILE='str_br1_4_ave.opt') ! output average hourly results for each day, Friant River outlet

      OPEN (UNIT=41,FILE='str_br1_1_max.opt')	! output maximum hourly results for each day, ML spillway
      OPEN (UNIT=42,FILE='str_br1_2_max.opt') ! output maximum hourly results for each day, F-K Canal
	OPEN (UNIT=43,FILE='str_br1_3_max.opt') ! output maximum hourly results for each day, Madera Canal
	OPEN (UNIT=44,FILE='str_br1_4_max.opt') ! output maximum hourly results for each day, Friant River outlet

      OPEN (UNIT=61,FILE='str_br1_1_min.opt')	! output minimum hourly results for each day, ML spillway
      OPEN (UNIT=62,FILE='str_br1_2_min.opt') ! output minimum hourly results for each day, F-K Canal
	OPEN (UNIT=63,FILE='str_br1_3_min.opt') ! output minimum hourly results for each day, Madera Canal
	OPEN (UNIT=64,FILE='str_br1_4_min.opt') ! output minimum hourly results for each day, Friant River outlet

      OPEN (UNIT=72,FILE='str_br1_ave.opt')   ! output average daily temperature results, average flow, ML
      OPEN (UNIT=74,FILE='str_wb1_vol.opt')   ! output average volume, ML

      OPEN (UNIT=82,FILE='str_br1_max.opt')   ! output daily, max hourly results, average flow, average volume, ML

      OPEN (UNIT=92,FILE='str_br1_min.opt')   ! output daily, min hourly results, average flow, average volume, ML

      OPEN (UNIT=101,FILE='output_1.opt')		! check
      OPEN (UNIT=102,FILE='output_2.opt')     ! check


      NOCT=0	
	NNOV=0
	NDEC=0
	NJAN=0
	NFEB=0
	NMAR=0
	NAPR=0
	NMAY=0
	NJUN=0
	NJUL=0
	NAUG=0
	NSEP=0

      SJDAY=21093. ! STARTING JDAY TO OUTPUT RESULTS (JDAY=21093 represents 10/1/1978)

      NYEAR=1978	! START FROM OCT. 1978, WATER YEAR 1979
     

      ! READ MILLERTON LAKE RESULTS

      NT_ML=0

	READ(2,*)
      READ(2,*)
15 	READ(2,*,END=20) TML1,TML2,TML3,TML4,TML5,TML6,TML7,TML8,TML9
     +                ,TML10,TML11,TML12,TML13

      IF ( TML1 .GE. SJDAY ) THEN

         NT_ML=NT_ML+1			! START FROM 1
	   TEMP_ML(1,NT_ML)=TML2	! ML SPILL TEMPERATURE
	   TEMP_ML(2,NT_ML)=TML3	! F-K CANAL TEMPERATURE
	   TEMP_ML(3,NT_ML)=TML4	! MADERA CANAL TEMPERATURE
	   TEMP_ML(4,NT_ML)=TML5	! FRIANT RIVER OUTLET TEMPERATURE
	   Q_ML(1,NT_ML)=TML6		! ML SPILL FLOW
	   Q_ML(2,NT_ML)=TML7		! F-K CANAL FLOW
	   Q_ML(3,NT_ML)=TML8		! MADERA CANAL FLOW
	   Q_ML(4,NT_ML)=TML9		! FRIANT RIVER OUTLET FLOW
	   E_ML(1,NT_ML)=TML10		! ML SPILL FLOW, ELEVATION
	   E_ML(2,NT_ML)=TML11		! F-K CANAL FLOW, ELEVATION
	   E_ML(3,NT_ML)=TML12		! MADERA CANAL FLOW, ELEVATION
	   E_ML(4,NT_ML)=TML13		! FRIANT RIVER OUTLET FLOW, ELEVATION

         GOTO 15

      ELSE

	   GOTO 15

      END IF

20    CONTINUE 

   

      ! READ MILLERTON LAKE VOLUME RESULTS

      NT_VML=0

      READ(4,*)
35 	READ(4,*,END=40) VML1,VML2,VML3,VML4

      IF ( VML1 .GE. SJDAY ) THEN
		
	   NT_VML=NT_VML+1			! START FROM 1
	   V_ML(1,NT_VML)=VML2		! ML TOTAL VOL
	   V_ML(2,NT_VML)=VML3		! ML VOL UNDER 1ST CRITERIA
	   V_ML(3,NT_VML)=VML4		! ML VOL UNDER 2ND CRITERIA

         GOTO 35

      ELSE 

	   GOTO 35

      END IF

40    CONTINUE 

c     check to see when the 1/24 day is off by one hr, the results show that JDAY=25758, it has 23 hours, this is a numerical issue

c      ttest=21093.004
c      ttest1=21093.004

c	do i=1,1000000

c	   if (ttest .lt. 30223.959 ) then

c		   ttest=ttest + 0.041667
c	       ttest1=ttest1 + (1./24.)
c		   write(103,'(2f10.3)') ttest,ttest1
		                
c		   if (int(ttest) .gt. int(ttest1) ) then
		   
c		       goto 999
			   
c		   end if	       

c         else
c             goto 999
c	   end if
         
c      end do

c999   continue

      write(*,*)  NT_ML,  NT_VML

c      stop


      DO I=1,3		 
         V_ML(I,NT_VML+1)=V_ML(I,NT_VML)         ! ENDING DAY HAS 23 HRS ONLY, NEED ONE MORE VALUE
	END DO

      DO I=1,4
         TEMP_ML(I,NT_ML+1)=TEMP_ML(I,NT_ML)	   ! ENDING DAY HAS 23 HRS ONLY, NEED ONE MORE VALUE 
         Q_ML(I,NT_ML+1)=Q_ML(I,NT_ML)
         E_ML(I,NT_ML+1)=E_ML(I,NT_ML)
	END DO




      DO I=1,3		
         V_ML(I,NT_VML+2)=V_ML(I,NT_VML)          ! total simulation periods have one hour off, NEED ONE MORE VALUE    
	END DO

      DO I=1,4
         TEMP_ML(I,NT_ML+2)=TEMP_ML(I,NT_ML)		! total simulation periods have one hour off, NEED ONE MORE VALUE 
         Q_ML(I,NT_ML+2)=Q_ML(I,NT_ML)
         E_ML(I,NT_ML+2)=E_ML(I,NT_ML)
	END DO



      !CHECK 

      WRITE(101,'("
     +   ML_SPILL_T    ML_F-KC_T    ML_MC_T    ML_RO_T
     +   ML_SPILL_Q    ML_F-KC_Q    ML_MC_Q    ML_RO_Q
     +   ML_SPILL_E    ML_F-KC_E    ML_MC_E    ML_RO_E   ")')

      DO I=1,223200

         WRITE(101,'(18F8.2)') 
     +                       TEMP_ML(1,I),TEMP_ML(2,I),
     +                       TEMP_ML(3,I),TEMP_ML(4,I),
     +                       Q_ML(1,I),Q_ML(2,I),
     +                       Q_ML(3,I),Q_ML(4,I),
     +                       E_ML(1,I),E_ML(2,I),
     +                       E_ML(3,I),E_ML(4,I)

	END DO


      ! DAILY RESULTS, CONSIDER (1) AVERAGE DAILY AND (2) MAXIMUM HOURLY FOR EACH DAY (3) MAXIMUM HOURLY FOR EACH DAY
	! (4) VARIATION (MAX MINUS MIN)

      NH=0
      TJDAY=21092	! 21092=9/30/1978, JDAY
               
	WRITE(72,'("    JDAY",<4>(4X,"AVET(C)"),<4>(2X,"AVEQ(cms)"),
     +     <4>(2X,"ELEV(m)"))')
                
	WRITE(74,'("    JDAY",(8X,"VOL_T(m3)"),(8X,"VOL_criteria(m3)"),
     +	      (8X,"VOL_criteria(m3)"))')	   
               
	WRITE(82,'("    JDAY",<4>(4X,"MAXT(C)"),<4>(2X,"AVEQ(cms)"))')   
                     
	WRITE(92,'("    JDAY",<4>(4X,"MINT(C)"),<4>(2X,"AVEQ(cms)"))')   
		 
	DO I=1,25

         NYEAR=NYEAR+1
         
	   DO J=1,31	   
	        
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0
			
		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			    
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH)) 
				   
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 				    

	           END IF
                  
              END DO
		  
	        NOCT=NOCT+1

		    TOCT_ML_1(1,NOCT)=SUM_ML1/24.
	        TOCT_ML_1(2,NOCT)=SUM_ML2/24. 
		    TOCT_ML_1(3,NOCT)=SUM_ML3/24.
	        TOCT_ML_1(4,NOCT)=SUM_ML4/24.
			
			TOCT_ML_2(1,NOCT)=TMAX_ML1
			TOCT_ML_2(2,NOCT)=TMAX_ML2
			TOCT_ML_2(3,NOCT)=TMAX_ML3
			TOCT_ML_2(4,NOCT)=TMAX_ML4    

			TOCT_ML_3(1,NOCT)=TMIN_ML1
			TOCT_ML_3(2,NOCT)=TMIN_ML2
			TOCT_ML_3(3,NOCT)=TMIN_ML3
			TOCT_ML_3(4,NOCT)=TMIN_ML4

		    QOCT_ML_1(1,NOCT)=SUM_QML1/24.
	        QOCT_ML_1(2,NOCT)=SUM_QML2/24. 
		    QOCT_ML_1(3,NOCT)=SUM_QML3/24.
	        QOCT_ML_1(4,NOCT)=SUM_QML4/24.

		    EOCT_ML_1(1,NOCT)=SUM_EML1/24.
	        EOCT_ML_1(2,NOCT)=SUM_EML2/24. 
		    EOCT_ML_1(3,NOCT)=SUM_EML3/24.
	        EOCT_ML_1(4,NOCT)=SUM_EML4/24.
			
	        TJDAY=TJDAY+1
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.

			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 

			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
          			               
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.

	   END DO
	   	            
         DO J=1,30
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			     
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))  
 
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NNOV=NNOV+1

		    TNOV_ML_1(1,NNOV)=SUM_ML1/24.
	        TNOV_ML_1(2,NNOV)=SUM_ML2/24. 
		    TNOV_ML_1(3,NNOV)=SUM_ML3/24.
	        TNOV_ML_1(4,NNOV)=SUM_ML4/24.
			
			TNOV_ML_2(1,NNOV)=TMAX_ML1
			TNOV_ML_2(2,NNOV)=TMAX_ML2
			TNOV_ML_2(3,NNOV)=TMAX_ML3
			TNOV_ML_2(4,NNOV)=TMAX_ML4 

			TNOV_ML_3(1,NNOV)=TMIN_ML1
			TNOV_ML_3(2,NNOV)=TMIN_ML2
			TNOV_ML_3(3,NNOV)=TMIN_ML3
			TNOV_ML_3(4,NNOV)=TMIN_ML4

		    QNOV_ML_1(1,NNOV)=SUM_QML1/24.
	        QNOV_ML_1(2,NNOV)=SUM_QML2/24. 
		    QNOV_ML_1(3,NNOV)=SUM_QML3/24.
	        QNOV_ML_1(4,NNOV)=SUM_QML4/24.

		    ENOV_ML_1(1,NNOV)=SUM_EML1/24.
	        ENOV_ML_1(2,NNOV)=SUM_EML2/24. 
		    ENOV_ML_1(3,NNOV)=SUM_EML3/24.
	        ENOV_ML_1(4,NNOV)=SUM_EML4/24.

	        TJDAY=TJDAY+1              
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.  			 
                
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
                
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
              
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

         DO J=1,31
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
			              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			     
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))  
 
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NDEC=NDEC+1

		    TDEC_ML_1(1,NDEC)=SUM_ML1/24.
	        TDEC_ML_1(2,NDEC)=SUM_ML2/24. 
		    TDEC_ML_1(3,NDEC)=SUM_ML3/24.
	        TDEC_ML_1(4,NDEC)=SUM_ML4/24.
			
			TDEC_ML_2(1,NDEC)=TMAX_ML1
			TDEC_ML_2(2,NDEC)=TMAX_ML2
			TDEC_ML_2(3,NDEC)=TMAX_ML3
			TDEC_ML_2(4,NDEC)=TMAX_ML4 

			TDEC_ML_3(1,NDEC)=TMIN_ML1
			TDEC_ML_3(2,NDEC)=TMIN_ML2
			TDEC_ML_3(3,NDEC)=TMIN_ML3
			TDEC_ML_3(4,NDEC)=TMIN_ML4

		    QDEC_ML_1(1,NDEC)=SUM_QML1/24.
	        QDEC_ML_1(2,NDEC)=SUM_QML2/24. 
		    QDEC_ML_1(3,NDEC)=SUM_QML3/24.
	        QDEC_ML_1(4,NDEC)=SUM_QML4/24.

		    EDEC_ML_1(1,NDEC)=SUM_EML1/24.
	        EDEC_ML_1(2,NDEC)=SUM_EML2/24. 
		    EDEC_ML_1(3,NDEC)=SUM_EML3/24.
	        EDEC_ML_1(4,NDEC)=SUM_EML4/24.

	        TJDAY=TJDAY+1

			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.  			 
                
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 

                
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
             
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

         DO J=1,31
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			     
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))   
  
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NJAN=NJAN+1

		    TJAN_ML_1(1,NJAN)=SUM_ML1/24.
	        TJAN_ML_1(2,NJAN)=SUM_ML2/24. 
		    TJAN_ML_1(3,NJAN)=SUM_ML3/24.
	        TJAN_ML_1(4,NJAN)=SUM_ML4/24.
			
			TJAN_ML_2(1,NJAN)=TMAX_ML1
			TJAN_ML_2(2,NJAN)=TMAX_ML2
			TJAN_ML_2(3,NJAN)=TMAX_ML3
			TJAN_ML_2(4,NJAN)=TMAX_ML4 

			TJAN_ML_3(1,NJAN)=TMIN_ML1
			TJAN_ML_3(2,NJAN)=TMIN_ML2
			TJAN_ML_3(3,NJAN)=TMIN_ML3
			TJAN_ML_3(4,NJAN)=TMIN_ML4

		    QJAN_ML_1(1,NJAN)=SUM_QML1/24.
	        QJAN_ML_1(2,NJAN)=SUM_QML2/24. 
		    QJAN_ML_1(3,NJAN)=SUM_QML3/24.
	        QJAN_ML_1(4,NJAN)=SUM_QML4/24.

		    EJAN_ML_1(1,NJAN)=SUM_EML1/24.
	        EJAN_ML_1(2,NJAN)=SUM_EML2/24. 
		    EJAN_ML_1(3,NJAN)=SUM_EML3/24.
	        EJAN_ML_1(4,NJAN)=SUM_EML4/24.

	        TJDAY=TJDAY+1               
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.  			 
                
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
                
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
              
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

         IF ( MOD(NYEAR,4) .EQ. 0 ) THEN
	   
	      DO J=1,29
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			     
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))   
 
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NFEB=NFEB+1

		    TFEB_ML_1(1,NFEB)=SUM_ML1/24.
	        TFEB_ML_1(2,NFEB)=SUM_ML2/24. 
		    TFEB_ML_1(3,NFEB)=SUM_ML3/24.
	        TFEB_ML_1(4,NFEB)=SUM_ML4/24.
			
			TFEB_ML_2(1,NFEB)=TMAX_ML1
			TFEB_ML_2(2,NFEB)=TMAX_ML2
			TFEB_ML_2(3,NFEB)=TMAX_ML3
			TFEB_ML_2(4,NFEB)=TMAX_ML4 

			TFEB_ML_3(1,NFEB)=TMIN_ML1
			TFEB_ML_3(2,NFEB)=TMIN_ML2
			TFEB_ML_3(3,NFEB)=TMIN_ML3
			TFEB_ML_3(4,NFEB)=TMIN_ML4

		    QFEB_ML_1(1,NFEB)=SUM_QML1/24.
	        QFEB_ML_1(2,NFEB)=SUM_QML2/24. 
		    QFEB_ML_1(3,NFEB)=SUM_QML3/24.
	        QFEB_ML_1(4,NFEB)=SUM_QML4/24.

		    EFEB_ML_1(1,NFEB)=SUM_EML1/24.
	        EFEB_ML_1(2,NFEB)=SUM_EML2/24. 
		    EFEB_ML_1(3,NFEB)=SUM_EML3/24.
	        EFEB_ML_1(4,NFEB)=SUM_EML4/24.

	        TJDAY=TJDAY+1               
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.		 
                
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
               
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
             
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	      END DO

         ELSE

	      DO J=1,28
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			     
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))  

			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 	

	           END IF
                  
              END DO
		  
	        NFEB=NFEB+1

		    TFEB_ML_1(1,NFEB)=SUM_ML1/24.
	        TFEB_ML_1(2,NFEB)=SUM_ML2/24. 
		    TFEB_ML_1(3,NFEB)=SUM_ML3/24.
	        TFEB_ML_1(4,NFEB)=SUM_ML4/24.
			
			TFEB_ML_2(1,NFEB)=TMAX_ML1
			TFEB_ML_2(2,NFEB)=TMAX_ML2
			TFEB_ML_2(3,NFEB)=TMAX_ML3
			TFEB_ML_2(4,NFEB)=TMAX_ML4 

			TFEB_ML_3(1,NFEB)=TMIN_ML1
			TFEB_ML_3(2,NFEB)=TMIN_ML2
			TFEB_ML_3(3,NFEB)=TMIN_ML3
			TFEB_ML_3(4,NFEB)=TMIN_ML4

		    QFEB_ML_1(1,NFEB)=SUM_QML1/24.
	        QFEB_ML_1(2,NFEB)=SUM_QML2/24. 
		    QFEB_ML_1(3,NFEB)=SUM_QML3/24.
	        QFEB_ML_1(4,NFEB)=SUM_QML4/24.

		    EFEB_ML_1(1,NFEB)=SUM_EML1/24.
	        EFEB_ML_1(2,NFEB)=SUM_EML2/24. 
		    EFEB_ML_1(3,NFEB)=SUM_EML3/24.
	        EFEB_ML_1(4,NFEB)=SUM_EML4/24.

	        TJDAY=TJDAY+1                
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.
                
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
               
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
              
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	      END DO

         END IF

         DO J=1,31
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			     
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))  
 
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NMAR=NMAR+1

		    TMAR_ML_1(1,NMAR)=SUM_ML1/24.
	        TMAR_ML_1(2,NMAR)=SUM_ML2/24. 
		    TMAR_ML_1(3,NMAR)=SUM_ML3/24.
	        TMAR_ML_1(4,NMAR)=SUM_ML4/24.
			
			TMAR_ML_2(1,NMAR)=TMAX_ML1
			TMAR_ML_2(2,NMAR)=TMAX_ML2
			TMAR_ML_2(3,NMAR)=TMAX_ML3
			TMAR_ML_2(4,NMAR)=TMAX_ML4 

			TMAR_ML_3(1,NMAR)=TMIN_ML1
			TMAR_ML_3(2,NMAR)=TMIN_ML2
			TMAR_ML_3(3,NMAR)=TMIN_ML3
			TMAR_ML_3(4,NMAR)=TMIN_ML4

		    QMAR_ML_1(1,NMAR)=SUM_QML1/24.
	        QMAR_ML_1(2,NMAR)=SUM_QML2/24. 
		    QMAR_ML_1(3,NMAR)=SUM_QML3/24.
	        QMAR_ML_1(4,NMAR)=SUM_QML4/24.

		    EMAR_ML_1(1,NMAR)=SUM_EML1/24.
	        EMAR_ML_1(2,NMAR)=SUM_EML2/24. 
		    EMAR_ML_1(3,NMAR)=SUM_EML3/24.
	        EMAR_ML_1(4,NMAR)=SUM_EML4/24.

	        TJDAY=TJDAY+1               
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.  			 
               
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
              
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
             
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

         DO J=1,30
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			      
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))  
  
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NAPR=NAPR+1

		    TAPR_ML_1(1,NAPR)=SUM_ML1/24.
	        TAPR_ML_1(2,NAPR)=SUM_ML2/24. 
		    TAPR_ML_1(3,NAPR)=SUM_ML3/24.
	        TAPR_ML_1(4,NAPR)=SUM_ML4/24.
			
			TAPR_ML_2(1,NAPR)=TMAX_ML1
			TAPR_ML_2(2,NAPR)=TMAX_ML2
			TAPR_ML_2(3,NAPR)=TMAX_ML3
			TAPR_ML_2(4,NAPR)=TMAX_ML4 

			TAPR_ML_3(1,NAPR)=TMIN_ML1
			TAPR_ML_3(2,NAPR)=TMIN_ML2
			TAPR_ML_3(3,NAPR)=TMIN_ML3
			TAPR_ML_3(4,NAPR)=TMIN_ML4

		    QAPR_ML_1(1,NAPR)=SUM_QML1/24.
	        QAPR_ML_1(2,NAPR)=SUM_QML2/24. 
		    QAPR_ML_1(3,NAPR)=SUM_QML3/24.
	        QAPR_ML_1(4,NAPR)=SUM_QML4/24.

		    EAPR_ML_1(1,NAPR)=SUM_EML1/24.
	        EAPR_ML_1(2,NAPR)=SUM_EML2/24. 
		    EAPR_ML_1(3,NAPR)=SUM_EML3/24.
	        EAPR_ML_1(4,NAPR)=SUM_EML4/24.

	        TJDAY=TJDAY+1               
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.  			 
               
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
                
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
            
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

         DO J=1,31
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			     
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))   
 
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NMAY=NMAY+1

		    TMAY_ML_1(1,NMAY)=SUM_ML1/24.
	        TMAY_ML_1(2,NMAY)=SUM_ML2/24. 
		    TMAY_ML_1(3,NMAY)=SUM_ML3/24.
	        TMAY_ML_1(4,NMAY)=SUM_ML4/24.
			
			TMAY_ML_2(1,NMAY)=TMAX_ML1
			TMAY_ML_2(2,NMAY)=TMAX_ML2
			TMAY_ML_2(3,NMAY)=TMAX_ML3
			TMAY_ML_2(4,NMAY)=TMAX_ML4 

			TMAY_ML_3(1,NMAY)=TMIN_ML1
			TMAY_ML_3(2,NMAY)=TMIN_ML2
			TMAY_ML_3(3,NMAY)=TMIN_ML3
			TMAY_ML_3(4,NMAY)=TMIN_ML4

		    QMAY_ML_1(1,NMAY)=SUM_QML1/24.
	        QMAY_ML_1(2,NMAY)=SUM_QML2/24. 
		    QMAY_ML_1(3,NMAY)=SUM_QML3/24.
	        QMAY_ML_1(4,NMAY)=SUM_QML4/24.

		    EMAY_ML_1(1,NMAY)=SUM_EML1/24.
	        EMAY_ML_1(2,NMAY)=SUM_EML2/24. 
		    EMAY_ML_1(3,NMAY)=SUM_EML3/24.
	        EMAY_ML_1(4,NMAY)=SUM_EML4/24.

	        TJDAY=TJDAY+1               
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.
                
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
               
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
              
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

         DO J=1,30
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			    
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))   

			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NJUN=NJUN+1

		    TJUN_ML_1(1,NJUN)=SUM_ML1/24.
	        TJUN_ML_1(2,NJUN)=SUM_ML2/24. 
		    TJUN_ML_1(3,NJUN)=SUM_ML3/24.
	        TJUN_ML_1(4,NJUN)=SUM_ML4/24.
			
			TJUN_ML_2(1,NJUN)=TMAX_ML1
			TJUN_ML_2(2,NJUN)=TMAX_ML2
			TJUN_ML_2(3,NJUN)=TMAX_ML3
			TJUN_ML_2(4,NJUN)=TMAX_ML4 

			TJUN_ML_3(1,NJUN)=TMIN_ML1
			TJUN_ML_3(2,NJUN)=TMIN_ML2
			TJUN_ML_3(3,NJUN)=TMIN_ML3
			TJUN_ML_3(4,NJUN)=TMIN_ML4

		    QJUN_ML_1(1,NJUN)=SUM_QML1/24.
	        QJUN_ML_1(2,NJUN)=SUM_QML2/24. 
		    QJUN_ML_1(3,NJUN)=SUM_QML3/24.
	        QJUN_ML_1(4,NJUN)=SUM_QML4/24.

		    EJUN_ML_1(1,NJUN)=SUM_EML1/24.
	        EJUN_ML_1(2,NJUN)=SUM_EML2/24. 
		    EJUN_ML_1(3,NJUN)=SUM_EML3/24.
	        EJUN_ML_1(4,NJUN)=SUM_EML4/24.

	        TJDAY=TJDAY+1              
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.  			 
               
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
               
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
             
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.

		  
	   END DO

         DO J=1,31
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			    
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))   

			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NJUL=NJUL+1

		    TJUL_ML_1(1,NJUL)=SUM_ML1/24.
	        TJUL_ML_1(2,NJUL)=SUM_ML2/24. 
		    TJUL_ML_1(3,NJUL)=SUM_ML3/24.
	        TJUL_ML_1(4,NJUL)=SUM_ML4/24.
			
			TJUL_ML_2(1,NJUL)=TMAX_ML1
			TJUL_ML_2(2,NJUL)=TMAX_ML2
			TJUL_ML_2(3,NJUL)=TMAX_ML3
			TJUL_ML_2(4,NJUL)=TMAX_ML4 

			TJUL_ML_3(1,NJUL)=TMIN_ML1
			TJUL_ML_3(2,NJUL)=TMIN_ML2
			TJUL_ML_3(3,NJUL)=TMIN_ML3
			TJUL_ML_3(4,NJUL)=TMIN_ML4

		    QJUL_ML_1(1,NJUL)=SUM_QML1/24.
	        QJUL_ML_1(2,NJUL)=SUM_QML2/24. 
		    QJUL_ML_1(3,NJUL)=SUM_QML3/24.
	        QJUL_ML_1(4,NJUL)=SUM_QML4/24.

		    EJUL_ML_1(1,NJUL)=SUM_EML1/24.
	        EJUL_ML_1(2,NJUL)=SUM_EML2/24. 
		    EJUL_ML_1(3,NJUL)=SUM_EML3/24.
	        EJUL_ML_1(4,NJUL)=SUM_EML4/24.

	        TJDAY=TJDAY+1               
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.			 
              
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
               
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
              
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

         DO J=1,31
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1
			   
			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			    
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))   
 
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NAUG=NAUG+1

		    TAUG_ML_1(1,NAUG)=SUM_ML1/24.
	        TAUG_ML_1(2,NAUG)=SUM_ML2/24. 
		    TAUG_ML_1(3,NAUG)=SUM_ML3/24.
	        TAUG_ML_1(4,NAUG)=SUM_ML4/24.
			
			TAUG_ML_2(1,NAUG)=TMAX_ML1
			TAUG_ML_2(2,NAUG)=TMAX_ML2
			TAUG_ML_2(3,NAUG)=TMAX_ML3
			TAUG_ML_2(4,NAUG)=TMAX_ML4 

			TAUG_ML_3(1,NAUG)=TMIN_ML1
			TAUG_ML_3(2,NAUG)=TMIN_ML2
			TAUG_ML_3(3,NAUG)=TMIN_ML3
			TAUG_ML_3(4,NAUG)=TMIN_ML4

		    QAUG_ML_1(1,NAUG)=SUM_QML1/24.
	        QAUG_ML_1(2,NAUG)=SUM_QML2/24. 
		    QAUG_ML_1(3,NAUG)=SUM_QML3/24.
	        QAUG_ML_1(4,NAUG)=SUM_QML4/24.

		    EAUG_ML_1(1,NAUG)=SUM_EML1/24.
	        EAUG_ML_1(2,NAUG)=SUM_EML2/24. 
		    EAUG_ML_1(3,NAUG)=SUM_EML3/24.
	        EAUG_ML_1(4,NAUG)=SUM_EML4/24.

	        TJDAY=TJDAY+1           
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.  			 
              
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
                
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
              
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

         DO J=1,30
	   
		    SUM_ML1=0.0
		    SUM_ML2=0.0
		    SUM_ML3=0.0
		    SUM_ML4=0.0

		    SUM_QML1=0.0
		    SUM_QML2=0.0
		    SUM_QML3=0.0
		    SUM_QML4=0.0

		    SUM_VML1=0.0
		    SUM_VML2=0.0
		    SUM_VML3=0.0

		    SUM_EML1=0.0
		    SUM_EML2=0.0
		    SUM_EML3=0.0
		    SUM_EML4=0.0
              
			DO K=1,24	   

			   NH=NH+1

			   SUM_ML1=SUM_ML1+TEMP_ML(1,NH)
			   SUM_ML2=SUM_ML2+TEMP_ML(2,NH)
			   SUM_ML3=SUM_ML3+TEMP_ML(3,NH)
			   SUM_ML4=SUM_ML4+TEMP_ML(4,NH)

			   SUM_QML1=SUM_QML1+Q_ML(1,NH)
			   SUM_QML2=SUM_QML2+Q_ML(2,NH)
			   SUM_QML3=SUM_QML3+Q_ML(3,NH)
			   SUM_QML4=SUM_QML4+Q_ML(4,NH)

			   SUM_VML1=SUM_VML1+V_ML(1,NH)
			   SUM_VML2=SUM_VML2+V_ML(2,NH)
			   SUM_VML3=SUM_VML3+V_ML(3,NH)

			   SUM_EML1=SUM_EML1+E_ML(1,NH)
			   SUM_EML2=SUM_EML2+E_ML(2,NH)
			   SUM_EML3=SUM_EML3+E_ML(3,NH)
			   SUM_EML4=SUM_EML4+E_ML(4,NH)

	           IF ( K.EQ. 1 ) THEN

                     TMAX_ML1=TEMP_ML(1,NH)
	               TMAX_ML2=TEMP_ML(2,NH)
                     TMAX_ML3=TEMP_ML(3,NH)
	               TMAX_ML4=TEMP_ML(4,NH)

                     TMIN_ML1=TEMP_ML(1,NH)
	               TMIN_ML2=TEMP_ML(2,NH)
                     TMIN_ML3=TEMP_ML(3,NH)
	               TMIN_ML4=TEMP_ML(4,NH)
                 
			   ELSE 
			     
			       TMAX_ML1=MAX(TMAX_ML1,TEMP_ML(1,NH))  
				   TMAX_ML2=MAX(TMAX_ML2,TEMP_ML(2,NH))  
				   TMAX_ML3=MAX(TMAX_ML3,TEMP_ML(3,NH))  
				   TMAX_ML4=MAX(TMAX_ML4,TEMP_ML(4,NH))  
 
			       TMIN_ML1=MIN(TMIN_ML1,TEMP_ML(1,NH))  
				   TMIN_ML2=MIN(TMIN_ML2,TEMP_ML(2,NH))  
				   TMIN_ML3=MIN(TMIN_ML3,TEMP_ML(3,NH))  
				   TMIN_ML4=MIN(TMIN_ML4,TEMP_ML(4,NH)) 

	           END IF
                  
              END DO
		  
	        NSEP=NSEP+1

		    TSEP_ML_1(1,NSEP)=SUM_ML1/24.
	        TSEP_ML_1(2,NSEP)=SUM_ML2/24. 
		    TSEP_ML_1(3,NSEP)=SUM_ML3/24.
	        TSEP_ML_1(4,NSEP)=SUM_ML4/24.
			
			TSEP_ML_2(1,NSEP)=TMAX_ML1
			TSEP_ML_2(2,NSEP)=TMAX_ML2
			TSEP_ML_2(3,NSEP)=TMAX_ML3
			TSEP_ML_2(4,NSEP)=TMAX_ML4 

			TSEP_ML_3(1,NSEP)=TMIN_ML1
			TSEP_ML_3(2,NSEP)=TMIN_ML2
			TSEP_ML_3(3,NSEP)=TMIN_ML3
			TSEP_ML_3(4,NSEP)=TMIN_ML4

		    QSEP_ML_1(1,NSEP)=SUM_QML1/24.
	        QSEP_ML_1(2,NSEP)=SUM_QML2/24. 
		    QSEP_ML_1(3,NSEP)=SUM_QML3/24.
	        QSEP_ML_1(4,NSEP)=SUM_QML4/24.

		    ESEP_ML_1(1,NSEP)=SUM_EML1/24.
	        ESEP_ML_1(2,NSEP)=SUM_EML2/24. 
		    ESEP_ML_1(3,NSEP)=SUM_EML3/24.
	        ESEP_ML_1(4,NSEP)=SUM_EML4/24.

	        TJDAY=TJDAY+1              
                     
			WRITE(72,'(13F10.2)') TJDAY,SUM_ML1/24.,SUM_ML2/24.,
     +                             SUM_ML3/24.,SUM_ML4/24.,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.,
     +                             SUM_EML1/24.,SUM_EML2/24.,
     +                             SUM_EML3/24.,SUM_EML4/24.	 
               
			WRITE(82,'(9F10.2)') TJDAY,TMAX_ML1,TMAX_ML2,TMAX_ML3,
     +			                 TMAX_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.     			 
                
			WRITE(92,'(9F10.2)') TJDAY,TMIN_ML1,TMIN_ML2,TMIN_ML3,
     +			                 TMIN_ML4,
     +                             SUM_QML1/24.,SUM_QML2/24.,
     +                             SUM_QML3/24.,SUM_QML4/24.
             
			WRITE(74,'(F10.2,3F16.2)') TJDAY,SUM_VML1/24.,SUM_VML2/24.
     +			                      ,SUM_VML3/24.
		  
	   END DO

	END DO


      WRITE(21,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	
	WRITE(22,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	 
	WRITE(23,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	
	WRITE(24,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'


      WRITE(41,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	
	WRITE(42,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	 
	WRITE(43,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	
	WRITE(44,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'


      WRITE(61,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	
	WRITE(62,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	 
	WRITE(63,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'
	
	WRITE(64,*) 'OCT   NOV   DEC   JAN   FEB   MAR   APR   MAY   JUN
     +   JUL   AUG   SEP'

      DO I=1,775

      ! MILLERTON LAKE
	   
	   WRITE(21,'(12F8.2)') TOCT_ML_1(1,I),TNOV_ML_1(1,I),
     + 	                    TDEC_ML_1(1,I),TJAN_ML_1(1,I),
     +                        TFEB_ML_1(1,I),TMAR_ML_1(1,I),
     +                        TAPR_ML_1(1,I),TMAY_ML_1(1,I),
     +                        TJUN_ML_1(1,I),TJUL_ML_1(1,I),
     +                        TAUG_ML_1(1,I),TSEP_ML_1(1,I)    

	   WRITE(22,'(12F8.2)') TOCT_ML_1(2,I),TNOV_ML_1(2,I),
     + 	                    TDEC_ML_1(2,I),TJAN_ML_1(2,I),
     +                        TFEB_ML_1(2,I),TMAR_ML_1(2,I),
     +                        TAPR_ML_1(2,I),TMAY_ML_1(2,I),
     +                        TJUN_ML_1(2,I),TJUL_ML_1(2,I),
     +                        TAUG_ML_1(2,I),TSEP_ML_1(2,I) 

	   WRITE(23,'(32F8.2)') TOCT_ML_1(3,I),TNOV_ML_1(3,I),
     + 	                    TDEC_ML_1(3,I),TJAN_ML_1(3,I),
     +                        TFEB_ML_1(3,I),TMAR_ML_1(3,I),
     +                        TAPR_ML_1(3,I),TMAY_ML_1(3,I),
     +                        TJUN_ML_1(3,I),TJUL_ML_1(3,I),
     +                        TAUG_ML_1(3,I),TSEP_ML_1(3,I)    

	   WRITE(24,'(12F8.2)') TOCT_ML_1(4,I),TNOV_ML_1(4,I),
     + 	                    TDEC_ML_1(4,I),TJAN_ML_1(4,I),
     +                        TFEB_ML_1(4,I),TMAR_ML_1(4,I),
     +                        TAPR_ML_1(4,I),TMAY_ML_1(4,I),
     +                        TJUN_ML_1(4,I),TJUL_ML_1(4,I),
     +                        TAUG_ML_1(4,I),TSEP_ML_1(4,I) 

	   WRITE(41,'(12F8.2)') TOCT_ML_2(1,I),TNOV_ML_2(1,I),
     + 	                    TDEC_ML_2(1,I),TJAN_ML_2(1,I),
     +                        TFEB_ML_2(1,I),TMAR_ML_2(1,I),
     +                        TAPR_ML_2(1,I),TMAY_ML_2(1,I),
     +                        TJUN_ML_2(1,I),TJUL_ML_2(1,I),
     +                        TAUG_ML_2(1,I),TSEP_ML_2(1,I)    

	   WRITE(42,'(12F8.2)') TOCT_ML_2(2,I),TNOV_ML_2(2,I),
     + 	                    TDEC_ML_2(2,I),TJAN_ML_2(2,I),
     +                        TFEB_ML_2(2,I),TMAR_ML_2(2,I),
     +                        TAPR_ML_2(2,I),TMAY_ML_2(2,I),
     +                        TJUN_ML_2(2,I),TJUL_ML_2(2,I),
     +                        TAUG_ML_2(2,I),TSEP_ML_2(2,I) 

	   WRITE(43,'(32F8.2)') TOCT_ML_2(3,I),TNOV_ML_2(3,I),
     + 	                    TDEC_ML_2(3,I),TJAN_ML_2(3,I),
     +                        TFEB_ML_2(3,I),TMAR_ML_2(3,I),
     +                        TAPR_ML_2(3,I),TMAY_ML_2(3,I),
     +                        TJUN_ML_2(3,I),TJUL_ML_2(3,I),
     +                        TAUG_ML_2(3,I),TSEP_ML_2(3,I)    

	   WRITE(44,'(12F8.2)') TOCT_ML_2(4,I),TNOV_ML_2(4,I),
     + 	                    TDEC_ML_2(4,I),TJAN_ML_2(4,I),
     +                        TFEB_ML_2(4,I),TMAR_ML_2(4,I),
     +                        TAPR_ML_2(4,I),TMAY_ML_2(4,I),
     +                        TJUN_ML_2(4,I),TJUL_ML_2(4,I),
     +                        TAUG_ML_2(4,I),TSEP_ML_2(4,I) 

	   WRITE(61,'(12F8.2)') TOCT_ML_3(1,I),TNOV_ML_3(1,I),
     + 	                    TDEC_ML_3(1,I),TJAN_ML_3(1,I),
     +                        TFEB_ML_3(1,I),TMAR_ML_3(1,I),
     +                        TAPR_ML_3(1,I),TMAY_ML_3(1,I),
     +                        TJUN_ML_3(1,I),TJUL_ML_3(1,I),
     +                        TAUG_ML_3(1,I),TSEP_ML_3(1,I)    

	   WRITE(62,'(12F8.2)') TOCT_ML_3(2,I),TNOV_ML_3(2,I),
     + 	                    TDEC_ML_3(2,I),TJAN_ML_3(2,I),
     +                        TFEB_ML_3(2,I),TMAR_ML_3(2,I),
     +                        TAPR_ML_3(2,I),TMAY_ML_3(2,I),
     +                        TJUN_ML_3(2,I),TJUL_ML_3(2,I),
     +                        TAUG_ML_3(2,I),TSEP_ML_3(2,I) 

	   WRITE(63,'(32F8.2)') TOCT_ML_3(3,I),TNOV_ML_3(3,I),
     + 	                    TDEC_ML_3(3,I),TJAN_ML_3(3,I),
     +                        TFEB_ML_3(3,I),TMAR_ML_3(3,I),
     +                        TAPR_ML_3(3,I),TMAY_ML_3(3,I),
     +                        TJUN_ML_3(3,I),TJUL_ML_3(3,I),
     +                        TAUG_ML_3(3,I),TSEP_ML_3(3,I)    

	   WRITE(64,'(12F8.2)') TOCT_ML_3(4,I),TNOV_ML_3(4,I),
     + 	                    TDEC_ML_3(4,I),TJAN_ML_3(4,I),
     +                        TFEB_ML_3(4,I),TMAR_ML_3(4,I),
     +                        TAPR_ML_3(4,I),TMAY_ML_3(4,I),
     +                        TJUN_ML_3(4,I),TJUL_ML_3(4,I),
     +                        TAUG_ML_3(4,I),TSEP_ML_3(4,I)

      END DO

      ! CHECK

      WRITE(102,'("
     +   ML_SPILL_T    ML_F-KC_T    ML_MC_T    ML_RO_T
     +   ML_SPILL_Q    ML_F-KC_Q    ML_MC_Q    ML_RO_Q
     +   ML_SPILL_E    ML_F-KC_E    ML_MC_E    ML_RO_E   ")')

      DO N=1,223200

         WRITE(102,'(18F8.2)') 
     +                       TEMP_ML(1,N),TEMP_ML(2,N),
     +                       TEMP_ML(3,N),TEMP_ML(4,N),
     +                       Q_ML(1,N),Q_ML(2,N),
     +                       Q_ML(3,N),Q_ML(4,N),
     +                       E_ML(1,N),E_ML(2,N),
     +                       E_ML(3,N),E_ML(4,N)

	END DO

      !===================================================================================================================
      ! CLOSE FILES

      CLOSE (UNIT = 2)

      CLOSE (UNIT = 4)

      CLOSE (UNIT = 21)
      CLOSE (UNIT = 22)
	CLOSE (UNIT = 23)
      CLOSE (UNIT = 24)

      CLOSE (UNIT = 41)
      CLOSE (UNIT = 42)
	CLOSE (UNIT = 43)
      CLOSE (UNIT = 44)

      CLOSE (UNIT = 61)
      CLOSE (UNIT = 62)
	CLOSE (UNIT = 63)
      CLOSE (UNIT = 64)

      CLOSE (UNIT = 72)
      CLOSE (UNIT = 74)

      CLOSE (UNIT = 82)

      CLOSE (UNIT = 92)

      CLOSE (UNIT = 101)
      CLOSE (UNIT = 102)

      END PROGRAM