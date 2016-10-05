PROGRAM writeHGplotData

   IMPLICIT NONE
   
   CHARACTER(LEN=1024) :: inType , vctFileN
   REAL(KIND=8) :: ScaleF , omega1
   INTEGER :: nOmega , nHG
   CHARACTER(LEN=1024) :: foldern , volin
   REAL(KIND=8) , DIMENSION(:) , ALLOCATABLE :: vctOmega

   REAL(KIND=8) , DIMENSIOn(:,:) , ALLOCATABLE :: matHG
   
   INTEGER :: IO
   INTEGER :: i1 , i2
   INTEGER , PARAMETER :: fid1 = 20 , fid2 = 21
   CHARACTER(LEN=1024) :: format_string , omegaString
   CHARACTER(LEN=1024) :: omfilen , outfilen
   COMPLEX(KIND=8) :: XVAL
   
   
   
   CALL readInput(inType,vctFileN,ScaleF,omega1,nOmega,nHG,foldern,volin,outfilen)
   
   IF     ( TRIM(ADJUSTL(inType)) == 'vector' ) THEN
      
      CALL readVctOmega(vctFileN,nOmega,vctOmega)
   
   ELSEIF ( TRIM(ADJUSTL(inType)) == 'equisp' ) THEN
      
      ALLOCATE(vctOmega(nOmega))
      
      vctOmega = (/ ( omega1 + (i1-1)*ScaleF , I1 = 1 , nOmega ) /)
      
   ENDIF
   
   !write(*,*) vctOmega
   ALLOCATE(matHG(SIZE(vctOmega),nHG))
   matHG = 0.0d0
   
   DO i1 =  1 , SIZE(vctOmega)
   
      format_string = "(A2,I2.2)"
      WRITE(omegaString,format_string) 'om' , i1
!      omfilen = TRIM(ADJUSTL(foldern)) // TRIM(ADJUSTL(omegaString)) // '/' // &
!                TRIM(ADJUSTL(volin)) // '/' // 'HarmonicGains_' // &
!                TRIM(ADJUSTL(volin)) // '_' // TRIM(ADJUSTL(omegaString)) // '.txt'
      omfilen = TRIM(foldern) // TRIM(omegaString) // 'HG.txt'

      OPEN(UNIT=FID2,FILE=omfilen)
      
      i2 = 0
      DO 
         IF ( I2 .LT. nHG) THEN
            READ(fid2,*,IOSTAT=io) xval  !matHG(i1,i2)
!            WRITE(*,*) io , i2
            IF ( io > 0) THEN
               WRITE(*,*) 'Check input.  Something was wrong'
               EXIT
            ELSE IF ( io < 0 ) THEN
               WRITE(*,*) 'EOF reached ! N. of read HGs = ' , i2
               EXIT
            ELSE
               I2 = I2 + 1
               matHG(i1,i2) = DBLE(xval)
            END IF
         ELSE
            EXIT
         END IF
      END DO
      CLOSE(FID2)
   END DO

!   OPEN(UNIT=FID1,FILE='./hgRe500.dat')
   OPEN(UNIT=FID1,FILE=TRIM(outfilen))
   WRITE(FID1,*) '    omega            HG(k)                   '
   WRITE(FID1,*) '---------------------------------------------------------------'
   DO I1 = 1 , SIZE(vctOmega)
!      DO I2 = 1 , nHG
         
         WRITE(FID1,*) vctOmega(i1) , matHG(i1,:)
         
!      END DO
   END DO
   
   CLOSE(FID1)
   
   
! ================================================================================   

   CONTAINS
   
   SUBROUTINE readInput(inType,vctFilen,ScaleF,omega1,nOmega,nHG,foldern,volin,outfilen)
   
   IMPLICIT NONE
   
   CHARACTER(LEN=1024) :: inType
   CHARACTER(LEN=1024) :: vctFilen
   REAL(KIND=8) :: ScaleF , omega1
   INTEGER :: nOmega , nHG
   CHARACTER(LEN=1024) :: foldern , volin , outfilen

   
   INTEGER , PARAMETER :: fid1 = 20
   
   OPEN(UNIT=fid1, FILE='./inputHgPlot.in')
   
   READ(fid1,*)
   READ(fid1,*) inType
   READ(fid1,*) vctFilen
   READ(fid1,*) ScaleF
   READ(fid1,*) omega1
   READ(fid1,*) nOmega
   READ(fid1,*)
   READ(fid1,*) nHG
   READ(fid1,*) foldern
   READ(fid1,*) volin
   READ(fid1,*) 
   READ(fid1,*) outfilen 

   
   CLOSE(fid1)
   
   END SUBROUTINE readInput

! ---------------------------------------------------------------------------------
  SUBROUTINE readVctOmega(vctFileN,nOmega,vctOmega)
   
   IMPLICIT NONE
   
   CHARACTER(LEN=1024) , INTENT(IN) :: vctFileN
   INTEGER      , INTENT(IN) :: nOmega
   REAL(KIND=8) , DIMENSION(:), ALLOCATABLE :: vctOmega
   
   REAL(KIND=8) :: omEl
   REAL(KIND=8) , DIMENSION(:), ALLOCATABLE :: vctOm
   
   INTEGER :: io , i1
   INTEGER , PARAMETER :: fid1 = 20
   
   ALLOCATE(vctOm(nOmega))
   vctOm = 0.0d0
   I1 = 0
   
   OPEN( UNIT=fid1, FILE=TRIM(ADJUSTL(vctFileN)) )
   
   DO 
      IF ( I1 .LT. nOmega) THEN
         READ(fid1,*,IOSTAT=io) omEl
         IF ( io > 0) THEN
            WRITE(*,*) 'Check input.  Something was wrong'
            EXIT
         ELSE IF ( io < 0 ) THEN
            WRITE(*,*) 'EOF reached ! N. of elements = ' , i1
            EXIT
         ELSE
            I1 = I1 + 1
            vctOm(i1) = omEl
         END IF
      ELSE
         EXIT
      END IF
   END DO
   
   CLOSE(fid1)
   
   ALLOCATE(vctOmega(i1))
   vctOmega = vctOm(1:i1)
   
   
   END SUBROUTINE readVctOmega

END PROGRAM writeHGplotData
