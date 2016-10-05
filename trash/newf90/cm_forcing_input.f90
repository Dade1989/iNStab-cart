MODULE cm_forcing_input

USE global_variables
USE sparse_matrix_operations
USE sparse_matrix_profiles
USE qv_sp
USE qc_sp_m
USE Dirichlet_Neumann
USE prep_mesh_p1p2_sp
USE start_sparse_kit
USE Gauss_points
USE vtk_plot

IMPLICIT NONE

CONTAINS

! ----------------------------------------------------------------------
  
  SUBROUTINE readForcingOrder ( a , bb  )
    
    IMPLICIT NONE
    
    INTEGER , DIMENSION(:)   , ALLOCATABLE :: a
    INTEGER , DIMENSION(:,:) , ALLOCATABLE :: bb
    INTEGER :: np
    INTEGER , PARAMETER :: nFor = 3
   
    np = 1
 
    ALLOCATE(a(nFor),bb(nFor,np))
    a(1)   = 0 ; bb(1,1) = 1 
    a(2)   = 2 ; bb(2,1) = 0 
    a(3)   = 1 ; bb(3,1) = 1 
    
    
  END SUBROUTINE readForcingOrder

! ----------------------------------------------------------------------

  SUBROUTINE switchForcing(x0,e0,str,ind,js_D,cm_Stiffness_Cmplx,f)
  
    REAL(KIND=8)    , DIMENSION(:)  , INTENT(IN) :: x0 , e0
    COMPLEX(KIND=8) , DIMENSION(:,:), INTENT(IN) :: str
    TYPE(dyn_int_line),DIMENSION(velCmpnnts) , INTENT(IN):: js_D
    COMPLEX(KIND=8) , DIMENSION(:,:), ALLOCATABLE :: str1 , str2
    INTEGER :: ind
    TYPE(CSR_MUMPS_Complex_Matrix)           :: cm_Stiffness_cmplx
    COMPLEX(KIND=8) , DIMENSION(:)  , ALLOCATABLE :: p1 , p2
    COMPLEX(KIND=8) , DIMENSION(:,:), ALLOCATABLE :: ff
    COMPLEX(KIND=8) , DIMENSION(:)  , ALLOCATABLE :: f
    COMPLEX(KIND=8) , DIMENSION(:,:), ALLOCATABLE :: fPlot
    
    COMPLEX(KIND=8) , DIMENSION(:)  , ALLOCATABLE :: x0_cmplx
    INTEGER :: i1
    IF ( ALLOCATED(f) ) THEN
      DEALLOCATE(f)
    ENDIF

    ALLOCATE(f(SIZE(x0,1)))
    f = 0.0d0
    
    IF ( ALLOCATED(x0_cmplx) ) THEN
      DEALLOCATE(x0_cmplx)
    ENDIF

    ALLOCATE( x0_cmplx(SIZE(x0)) )
    x0_cmplx = x0
    
    SELECT CASE(ind)
    CASE(1)
! p(1) < Dv . DU >         + b.c.   MOVED TO a0bbForcing(...) below
!      CALL zAtimx( f , cm_Stiffness_Cmplx % e , &
!                       cm_Stiffness_Cmplx % j , &
!                       cm_Stiffness_Cmplx % i , &
!                       x0_cmplx )
!  
    CASE(2)
      ! < v . [(u.D)u] >         + b.c.
      ALLOCATE(str1(np,k_d),str2(np,k_d)) 
      CALL extract_Cmplx (str(:,1),str1,p1) 
      CALL extract_Cmplx (str(:,2),str2,p2) 
      p1 = 0.0d0
  

  
      ALLOCATE(ff(k_d,np))
      ff = 0.0d0
      CALL qv_001_sp_Cmplx ( mm,jj,str1,str2, ff )
      CALL collect_cmplx (-ff,p1,f)
      ! DEALLOCATE(ff)
      ! Impose boundary conditions
      CALL Dirichlet_homo_c_Cmplx(np,js_D,f)

!   ! Check ----
!         CALL vtk_plot_eigenvectors (rr, jj,  &
!                DBLE(str(:,1:1)), &
!                trim(p_in%plot_directory)// &
!                'str1.vtk')
!   ! ----------
!   ! Check ----
!         CALL vtk_plot_eigenvectors (rr, jj,  &
!                DBLE(str(:,2:2)), &
!                trim(p_in%plot_directory)// &
!                'str2.vtk')
!   ! ----------
!   ! Check ----
!         IF (ALLOCATED(fPlot)) THEN ;  DEALLOCATE(fPlot) ; ENDIF
!         ALLOCATE(fPlot(SIZE(f),1))
!         fPlot(:,1) = f
!         CALL vtk_plot_eigenvectors (rr, jj,  &
!                DBLE(fPlot(:,1:1)), &
!                trim(p_in%plot_directory)// &
!                'gDg.vtk')
!   ! ----------

    CASE(3)
      ! p(1) < Dv . Du >         + b.c.
      CALL zAtimx( f , cm_Stiffness_Cmplx % e , &
                       cm_Stiffness_Cmplx % j , &
                       cm_Stiffness_Cmplx % i , &
                       str(:,1) )
      CALL Dirichlet_homo_c_Cmplx(np,js_D,f)

!   ! Check ----
!   IF (ALLOCATED(fPlot)) THEN ;  DEALLOCATE(fPlot) ; ENDIF
!   ALLOCATE(fPlot(SIZE(f),1))
!   fPlot(:,1) = f
!   OPEN(UNIT=20,FILE="./plots/f_check.txt")
!   DO i1 = 1 , np
!     WRITE(20,*) f(i1) , f(i1+np)
!   ENDDO
!   CLOSE(20)
!         CALL vtk_plot_eigenvectors (rr, jj,  &
!                DBLE(fPlot(:,1:1)), &
!                trim(p_in%plot_directory)// &
!                'Ku.vtk')
!   STOP
!   ! ----------
    CASE DEFAULT
      WRITE(*,*) "Error in switchForcing() : Wrong input parameter"
      WRITE(*,*)
      STOP
    ENDSELECT
  
  END SUBROUTINE switchForcing

! ----------------------------------------------------------------------

SUBROUTINE a0bbForcing ( x0,e0,cm_Stiffness_Cmplx,bb, f )

    REAL(KIND=8)    , DIMENSION(:)  , INTENT(IN) :: x0 , e0
    TYPE(CSR_MUMPS_Complex_Matrix)  , INTENT(IN) :: cm_Stiffness_cmplx
    INTEGER , DIMENSION(:)          , INTENT(IN) :: bb
    COMPLEX(KIND=8) , DIMENSION(:)  , ALLOCATABLE :: f
    COMPLEX(KIND=8) , DIMENSION(:,:), ALLOCATABLE :: fPlot
    
    COMPLEX(KIND=8) , DIMENSION(:)  , ALLOCATABLE :: x0_cmplx

    INTEGER :: i1 

    IF ( ALLOCATED(f) ) THEN
      DEALLOCATE(f)
    ENDIF

    ALLOCATE(f(SIZE(x0,1)))
    f = 0.0d0
    
    IF ( ALLOCATED(x0_cmplx) ) THEN
      DEALLOCATE(x0_cmplx)
    ENDIF

    ALLOCATE( x0_cmplx(SIZE(x0)) )
    x0_cmplx = x0


    IF ( ANY(bb.EQ.(/1/)) ) THEN
!      WRITE(*,*) " I entered a0bbForcing () "
      ! p(1) < Dv . DU >         + b.c.
      CALL zAtimx( f , cm_Stiffness_Cmplx % e , &
                       cm_Stiffness_Cmplx % j , &
                       cm_Stiffness_Cmplx % i , &
                       x0_cmplx )
!      WRITE(*,*) "max(abs(f01)) = " , MAXVAL(ABS(f))
      CALL Dirichlet_homo_c_Cmplx(np,js_D,f)
    ENDIF

!   ! Check ----
!   WRITE(*,*) "SIZE(jj,2) = " , SIZE(jj,2)
!   OPEN(UNIT=20,FILE="./plots/f_el.txt")
!   DO i1 = 1 , SIZE(jj,2)
!    WRITE(20,*) f(jj(:,i1)) 
!   ENDDO
!   CLOSE(20)
!   STOP
!   
!   IF (ALLOCATED(fPlot)) THEN ;  DEALLOCATE(fPlot) ; ENDIF
!   ALLOCATE(fPlot(SIZE(f),1))
!   fPlot(:,1) = f
!   OPEN(UNIT=20,FILE="./plots/f_check.txt")
!   DO i1 = 1 , np
!     WRITE(20,*) f(i1) , f(i1+np)
!   ENDDO
!   CLOSE(20)
!         CALL vtk_plot_eigenvectors (rr, jj,  &
!                DBLE(fPlot(:,1:1)), &
!                trim(p_in%plot_directory)// &
!                'Ku001.vtk')
!   STOP
!   ! ----------

END SUBROUTINE a0bbForcing

! ----------------------------------------------------------------------

END MODULE cm_forcing_input
