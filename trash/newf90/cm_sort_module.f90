MODULE cm_sort_module

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module containing the SUBROUTINEs to sort vector (only INTEGER, to now)
!
!
! 2016-03-15
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! SUBROUTINES:
!  sort            ( order , vSort , iSort )
!  flip            ( v )
!  ascendSort      ( vSort , iSort )
! 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IMPLICIT NONE

CONTAINS

! ---------------------------------------------------------------------------

SUBROUTINE sort(order,vSort,iSort)
! only for integer (to now)
! 
! Input :
!  vec : vector of integer
!  order : 'ascending' , 
! Output :
!  

IMPLICIT NONE

CHARACTER(*) , INTENT(IN) :: order
INTEGER , DIMENSION(:) , INTENT(INOUT) :: vSort
INTEGER , DIMENSION(:) , ALLOCATABLE :: iSort
INTEGER , DIMENSION(:) , ALLOCATABLE :: vec

INTEGER :: iOrd
INTEGER , DIMENSION(:) , ALLOCATABLE :: a_d
INTEGER :: n_a_d

ALLOCATE(vec(SIZE(vSort)))
vec = vSort

! Check on the inputs
IF ( order == 'ascend' ) THEN
 iOrd = 0;
ELSEIF ( order == 'descend' ) THEN
 iOrd = 1;
ELSE
 WRITE(*,*) " "
 WRITE(*,*) " error in sort (...,order,...,...)      "
 WRITE(*,*) " --> order must be 'ascend' or 'descend'"
 WRITE(*,*) " "
 RETURN
ENDIF

! order = 'ascending'
 CALL ascendSort(vSort,iSort)
IF ( order == 'descend' ) THEN
  CALL flip(vSort)
  CALL flip(iSort)
END IF 


END SUBROUTINE sort

! ---------------------------------------------------------------------------

SUBROUTINE flip (v)
! only for integer (to now)

IMPLICIT NONE

INTEGER , DIMENSION(:) , INTENT(INOUT) :: v
INTEGER , DIMENSION(:) , ALLOCATABLE :: v0
INTEGER :: n
INTEGER :: i1

n = SIZE(v,1)
ALLOCATE(v0(n))
v0 = v

DO i1 = 1 , n
  v(i1) = v0(n+1-i1) 
END DO

END SUBROUTINE flip

! ---------------------------------------------------------------------------

SUBROUTINE ascendSort ( vSort , iSort )
   
   IMPLICIT NONE
   
   INTEGER , DIMENSION(:) , INTENT(INOUT) :: vSort
   INTEGER , DIMENSION(:) , ALLOCATABLE   :: iSort
   
   INTEGER ,  DIMENSION(:) , ALLOCATABLE :: A
   INTEGER :: M , MINOLD , EL
   INTEGER :: MINV , KMIN
   INTEGER :: I1 , i2
   
   ALLOCATE( ISORT(SIZE(VSORT)) , A(SIZE(VSORT)) )
   ISORT = 0
   A     = VSORT
   
   M    = MAXVAL(A) + 1
   MINV = MINVAL(A) - 1
   EL = 0
   DO I1 = 1 , SIZE(A)
      
      MINOLD = MINV
      MINV = A(1)
      KMIN = 1
      DO I2 = 1 , SIZE(A)
         
         IF ( A(I2) .LE. MINV ) THEN
            MINV = A(I2)
            KMIN = I2
            CYCLE
         END IF
         
      END DO
      IF ( MINV .NE. MINOLD ) THEN
         EL = EL + 1
      END IF
      VSORT(I1) = A(KMIN)
      ISORT(I1) = KMIN
      A(KMIN)   = M
   END DO


END SUBROUTINE ascendSort


! ----------------------------------------------------------------------------

END MODULE cm_sort_module
