MODULE cm_time_evolution

 USE global_variables

 IMPLICIT NONE
 
 CONTAINS

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 SUBROUTINE evolution_on_centre_manifold()
 
 IMPLICIT NONE
 
   INTEGER :: nc , npa
   INTEGER         , DIMENSION(:,:) , ALLOCATABLE :: matIndTot
   COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: GG
   COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: a0
   COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: at
   
   REAL(KIND=8) :: eRe

   COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: cFt
   
   ! Read matIndTot and GG coefficients
   CALL readGG ( nc , npa , matIndTot , GG)

   ! Bifurcation parameter : e = (Re-Re_c)/(Re*Re_c)
   eRe = ( p_in%cm_Re - Re) / (p_in%cm_Re*Re)
!  ! Check ----
!  WRITE(*,*) p_in%cm_Re
!  WRITE(*,*) Re
!  WRITE(*,*) eRe
!  ! ----------
   ! Initial condition :  compute it! ( cyl: It must be c.c. )
   ! ...
   ALLOCATE(a0(nc))
!  Order 3
!   a0(1) = CMPLX(1.76892d0,0.0d0,KIND=8) 
!   a0(2) = CMPLX(1.76892d0,0.0d0,KIND=8) 
!  Order 4
   a0(1) = CMPLX(1.8605d0,0.0d0,KIND=8) 
   a0(2) = CMPLX(1.8605d0,0.0d0,KIND=8) 

   a0(1) = CMPLX(1.0d0,0.0d0,KIND=8) 
   a0(2) = CMPLX(1.0d0,0.0d0,KIND=8) 

!   CALL complexRK4 ( matIndTot , GG , (/ eRe /) , a0 , at )
 
!   CALL force_coefficients_evolution ( matIndTot , at , (/ eRe /) , cFt )
   WRITE(*,*) " eRe = " , eRe
   CALL limit_cycle_solution ( matIndTot , eRe , 4 )
 
 END SUBROUTINE evolution_on_centre_manifold
 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 SUBROUTINE complexRK4 ( matIndTot , GG , eRe , a0 , at )
 
   IMPLICIT NONE
   
   INTEGER         , DIMENSION(:,:) , INTENT(IN) :: matIndTot
   COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: GG
   REAL(KIND=8)    , DIMENSION(:)   , INTENT(IN) :: eRe
   COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN) :: a0
   
   COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: at
   
   REAL(KIND=8)    , DIMENSION(:)   , ALLOCATABLE :: aTab , bTab
   COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: k1 , k2 , k3 , k4
   COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: a1
   REAL(KIND=8) :: t0 , tF , dt
   INTEGER      :: ndt
   
   INTEGER :: i1
   
   ! Time parameters
   t0 = p_in%cm_tInit
   tF = p_in%cm_tEnd
   dt = p_in%cm_dt
   
   ndt = CEILING( (tF-t0)/dt ) 
   
   ! Runge Kutta coefficients
   ALLOCATE(at(SIZE(a0),ndt+1))
   ALLOCATE(aTab(3),bTab(4))
   aTab(1) = 1.0d0 / 2.0d0;aTab(2) = 1.0d0 / 2.0d0;aTab(3) = 1.0d0
   
   bTab(1) = 1.0d0 / 6.0d0;bTab(2) = 1.0d0 / 3.0d0
   bTab(3) = 1.0d0 / 3.0d0;bTab(4) = 1.0d0 / 6.0d0
   
   ! Initial conditions
   at(:,1) = a0
   
   ALLOCATE(k1(SIZE(a0)))
   ALLOCATE(k2(SIZE(a0)))
   ALLOCATE(k3(SIZE(a0)))
   ALLOCATE(k4(SIZE(a0)))
   ALLOCATE(a1(SIZE(a0)))
   
   
   ! Time integration
   DO i1 = 1 , ndt
     
     CALL forcingG( matIndTot , GG , at(:,i1)               , eRe , k1 )
     CALL forcingG( matIndTot , GG , at(:,i1)+dt*aTab(1)*k1 , eRe , k2 )
     CALL forcingG( matIndTot , GG , at(:,i1)+dt*aTab(2)*k2 , eRe , k3 )
     CALL forcingG( matIndTot , GG , at(:,i1)+dt*aTab(3)*k3 , eRe , k4 )
     
     at(:,i1+1) = at(:,i1) + dt * ( bTab(1) * k1 + bTab(2) * k2 + & 
                                    bTab(3) * k3 + bTab(4) * k4   )

   ENDDO 
   
   
   OPEN(UNIT=24,FILE='./plots/aa.txt')
   DO i1 = 1 , SIZE(at,2)
     WRITE(24,*) at(:,i1)
   ENDDO
   CLOSE(24)

   OPEN(UNIT=24,FILE='./plots/aaReal.txt')
   DO i1 = 1 , SIZE(at,2)
     WRITE(24,*) DBLE(at(1,i1)) , AIMAG(at(1,i1)) , &
                 DBLE(at(2,i1)) , AIMAG(at(2,i1)) , &
                 SQRT( DBLE(at(1,i1))* DBLE(at(2,i1))- &
                      AIMAG(at(1,i1))*AIMAG(at(2,i1)) )
   ENDDO
   CLOSE(24)
 
 END SUBROUTINE complexRK4
 
 
 ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 SUBROUTINE forcingG ( matIndTot , GG , a , e , f )
 
   IMPLICIT NONE
   
   INTEGER         , DIMENSION(:,:) , INTENT(IN) :: matIndTot
   COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: GG
   COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN) :: a
   REAL(KIND=8)    , DIMENSION(:)   , INTENT(IN) :: e
   
   COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: f
   
   INTEGER :: i1
   
   IF ( ALLOCATED(f) ) DEALLOCATE(f)
   ALLOCATE(f(SIZE(a)))
   
   f = CMPLX(0.0d0,0.0d0,KIND=8)

!WRITE(*,*) " a = " , a
!WRITE(*,*) " e = " , e
   
   DO i1 = 1 ,  SIZE(matIndTot,1) ! 19
!WRITE(*,*) a
!WRITE(*,*) matIndTot(i1,1:SIZE(a))
     f = f + GG(:,i1) * PRODUCT( a**matIndTot(i1,1:SIZE(a)) ) * & 
                        PRODUCT( CMPLX( e,0.0d0,KIND=8) &
                         **matIndTot(i1,SIZE(a)+1:SIZE(a)+SIZE(e)) ) 

   ENDDO
   
 
 END SUBROUTINE forcingG

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 SUBROUTINE readGG( nc , npa , matIndTot , GG)
 
   IMPLICIT NONE
   
   INTEGER :: nc , npa , nTerms
   INTEGER , DIMENSION(:,:) , ALLOCATABLE :: matIndTot
   COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: GG
   
   INTEGER :: i1 
   
   OPEN(UNIT=23,FILE="./plots/GG.txt")
   READ(23,*)
   READ(23,*) nc , npa , nTerms
   READ(23,*)
   
   ALLOCATE(matIndTot(nTerms,nc+npa))
   ALLOCATE(GG(nc,nTerms))
   DO i1 = 1 , nTerms
WRITE(*,*)  i1 
     READ(23,*) matIndTot(i1,:) , GG(:,i1) 
   
   ENDDO
   CLOSE(23)
 
 END SUBROUTINE readGG
 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE  force_coefficients_evolution ( matIndTot , a , e , cFt )
 
   IMPLICIT NONE
   
   INTEGER         , DIMENSION(:,:) , INTENT(IN) :: matIndTot
   COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: a
   REAL   (KIND=8) , DIMENSION(:)   , INTENT(IN) :: e
   COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: cFt
   
   COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: cFcm
   REAL(KIND=8) :: rcFcm1 , icFcm1 , rcFcm2 , icFcm2
   COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: f

   REAL(KIND=8) :: tol   
   INTEGER :: i1 ,i2 , dummy , iplot 

   ALLOCATE(cFcm(2,size(matIndTot,1)))
   OPEN(UNIT=21,FILE='./plots/cmForceCoeff.txt')
   READ(21,*)
   READ(21,*)
   READ(21,*)
   DO i1 = 1 , size(matIndTot,1)
     READ(21,*) dummy , rcFcm1 , icFcm1, rcFcm2 , icFcm2
     cFcm(1,i1) = CMPLX(rcFcm1,icFcm1,KIND=8)
     cFcm(2,i1) = CMPLX(rcFcm2,icFcm2,KIND=8)
   ENDDO
   CLOSE(21)
! Check ----
WRITE(*,*)
DO i1 = 1 , SIZE(cFcm,2)
  WRITE(*,*) i1 , cFcm(:,i1) , matIndTot(i1,:)
ENDDO
WRITE(*,*)
!-----------

! for coarse meshes ----
! CALL cleanCf ( cFcm )
! ----------------------
   
   IF ( ALLOCATED(cFt) ) DEALLOCATE(cFt)
   ALLOCATE(cFt(2,SIZE(a,2)))

   cft = CMPLX(0.0d0,0.0d0,KIND=8)

WRITE(*,*) " SIZE(matIndTot,1) = " , SIZE(matIndTot,1)
WRITE(*,*) " eRe               = " , e              
   DO i2 = 1 , SIZE(a,2)
     DO i1 = 1 , SIZE(matIndTot,1) !19
 
      cft(1,i2) = cft(1,i2) &
                + cFcm(1,i1) * ( a(1,i2) ** matIndTot(i1,1) * & 
                                 a(2,i2) ** matIndTot(i1,2) * & 
                                 e(1)    ** matIndTot(i1,3)  )

      cft(2,i2) = cft(2,i2) &
                + cFcm(2,i1) * ( a(1,i2) ** matIndTot(i1,1) * & 
                                 a(2,i2) ** matIndTot(i1,2) * & 
                                 e(1)    ** matIndTot(i1,3)  )
  ! PRODUCT( a(:,i2)**matIndTot(i1,1:SIZE(a,1))) * & 
  !                             PRODUCT( e  &
  !                           **matIndTot(i1,SIZE(a,1)+1:SIZE(a,1)+SIZE(e))) 
 
     ENDDO
   ENDDO


   iplot = 1
   IF ( SIZE(a,2) .GT. 100001 ) THEN
      iplot = SIZE(a,2) / 100000
   ENDIF
   WRITE(*,*) "iplot = " , iplot

   tol = 1e-15
   OPEN(UNIT=22,FILE='./plots/cmFt.txt')
   WRITE(22,*) "#          cD                       cL        "
   IF ( ( MAXVAL(ABS(AIMAG(cft(1,:)))) .LT. tol ) .AND. &
        ( MAXVAL(ABS(AIMAG(cft(2,:)))) .LT. tol ) )    THEN
     WRITE(*,*)
     WRITE(*,*) " Force coefficients are real within tolerance"
     WRITE(*,*) " tol = " , tol 
     WRITE(*,*)

     DO i1 = 1 , SIZE(cft,2)
       IF ( MOD(i1,iplot) .EQ. 0 ) THEN
         WRITE(22,*) DBLE(cft(1,i1)) , DBLE(cft(2,i1))
       ENDIF
     ENDDO


   ELSE
     WRITE(*,*)
     WRITE(*,*) " Error : force coefficients are complex!    "
     WRITE(*,*) "         Something may be wrong.            "
     WRITE(*,*) "         Complex output in ./plots/cmFt.txt "
     WRITE(*,*) " MAXVAL(ABS(AIMAG(CD))) = " , MAXVAL(ABS(AIMAG(CFT(1,:))))
     WRITE(*,*) " MAXVAL(ABS(AIMAG(CL))) = " , MAXVAL(ABS(AIMAG(CFT(2,:))))
     WRITE(*,*) " tol = " , tol 
     WRITE(*,*)

     DO i1 = 1 , SIZE(cft,2)
       IF ( MOD(i1,iplot) .EQ. 0 ) THEN
         WRITE(22,*) cft(:,i1)
       ENDIF
     ENDDO

   ENDIF
   CLOSE(22)

   
END SUBROUTINE force_coefficients_evolution

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE limit_cycle_solution ( matIndTot , eRe , order )

IMPLICIT NONE

INTEGER , DIMENSION(:,:) , INTENT(IN) :: matIndTot
INTEGER ,                  INTENT(IN) :: order
REAL(KIND=8) :: eRe


INTEGER :: nc , npa
INTEGER , DIMENSION(:,:) , ALLOCATABLE :: matIndTot1
COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: GG
COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: aLc , cfLc 
REAL(KIND=8) :: rLc , oLc , period , ti , dt

COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: cFt

INTEGER :: i1 , ndt


nc  = 2
npa = 1

! Read GG
CALL readGG( nc , npa , matIndTot1 , GG)

! Find the radius and the pulastion on the LC
IF ( order .EQ. 3 ) THEN
  rLc = SQRT( - ( DBLE(GG(1,7)) *eRe + eRe**2 * DBLE(GG(1,17)) ) / &
                ( DBLE(GG(1,11)) ) )
  oLc = AIMAG(GG(1,1 ))                + &
        AIMAG(GG(1,7 )) * eRe          + &
        AIMAG(GG(1,11)) * rLc**2       + &
        AIMAG(GG(1,17)) * eRe**2
ELSEIF ( order .EQ. 4 ) THEN
  rLc = SQRT(- ( DBLE(GG(1,7))*eRe + eRe**2 * DBLE(GG(1,17)) + &
                 DBLE(GG(1,32))* eRe**3 ) / &
               ( DBLE(GG(1,11)) + DBLE(GG(1,26))*eRe ) )
  oLc = AIMAG(GG(1,1 ))                + &
        AIMAG(GG(1,7 )) * eRe          + &
        AIMAG(GG(1,11)) * rLc**2       + &
        AIMAG(GG(1,17)) * eRe**2       + &
        AIMAG(GG(1,26)) * rLc**2 * eRe + &
        AIMAG(GG(1,32)) * eRe**3
ELSE 
  WRITE(*,*) " Error in limit_cycle_solution () "
  WRITE(*,*) " for order > 4 , not implemented to now!"
  STOP
ENDIF

period = 2.0d0 * 4.0d0 * ATAN(1.0d0)  / oLc
WRITE(*,*)
WRITE(*,*) "eRe       = " , eRe
WRITE(*,*) "r_LC      = " , rLC
WRITE(*,*) "omega_LC  = " , oLC
WRITE(*,*) "period_LC = " , period
WRITE(*,*)

! Write the periodic solution on the LC
ndt = 100
dt  = period/DBLE(ndt)
WRITE(*,*) "dt = " , dt 
ALLOCATE(aLC(2,ndt+1))
DO i1 = 0 , ndt
  ti = i1 * dt 
  aLc(1,i1+1) = rLc * CMPLX( COS(oLc*ti) , SIN(oLc*ti) , KIND=8 )
  aLc(2,i1+1) = rLc * CMPLX( COS(oLc*ti) ,-SIN(oLc*ti) , KIND=8 )
ENDDO

   OPEN(UNIT=24,FILE='./plots/aaLCreal.txt')
   DO i1 = 1 , SIZE(aLC,2)
     WRITE(24,*) DBLE(aLC(1,i1)) , AIMAG(aLc(1,i1)) , &
                 DBLE(aLC(2,i1)) , AIMAG(aLc(2,i1))    
   ENDDO
   CLOSE(24)

CALL force_coefficients_evolution ( matIndTot , aLC , (/eRe/) , cFt )


END SUBROUTINE limit_cycle_solution

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE cleanCf ( c )

IMPLICIT NONE

COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(INOUT) :: c

WRITE(*,*) SIZE(c,1) , SIZE(c,2)

c(1,1 ) = 0.0d0
c(1,2 ) = 0.0d0
c(2,3 ) = 0.0d0
c(2,4 ) = 0.0d0
c(2,5 ) = 0.0d0
c(2,6 ) = 0.0d0
c(1,7 ) = 0.0d0
c(1,8 ) = 0.0d0
c(2,9 ) = 0.0d0
c(1,10) = 0.0d0
c(1,11) = 0.0d0
c(1,12) = 0.0d0
c(1,13) = 0.0d0
c(2,14) = 0.0d0
c(2,15) = 0.0d0
c(2,16) = 0.0d0
c(1,17) = 0.0d0
c(1,18) = 0.0d0
c(2,19) = 0.0d0
c(2,20) = 0.0d0
c(2,21) = 0.0d0
c(2,22) = 0.0d0
c(2,23) = 0.0d0
c(2,24) = 0.0d0
c(1,25) = 0.0d0
c(1,26) = 0.0d0
c(1,27) = 0.0d0
c(1,28) = 0.0d0
c(2,29) = 0.0d0
c(2,30) = 0.0d0
c(2,31) = 0.0d0
c(1,32) = 0.0d0
c(1,33) = 0.0d0 
c(2,34) = 0.0d0

END SUBROUTINE cleanCf

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
END MODULE cm_time_evolution
