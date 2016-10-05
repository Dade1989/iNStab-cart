MODULE hg_compute_hg_sens

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! HARMONIC GAIN SENSITIVTY
!  only code 2d TO NOW: wavenumber = 0 in z direction
!
! 
! 2016-05-12
! 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SUBROUTINES : 
!  compute_hg_adjoint           ( )
!  compute_hgSens_baseFlow      ( )
!  compute_hgSens_steadyControl ( )
!  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  USE par_solve_mumps
  USE sparse_matrix_profiles
  USE sparse_matrix_operations
  USE prep_mesh_p1p2_sp
  USE global_variables
  USE qc_sp_m
  USE vtk_plot
  USE fem_miscellaneous

IMPLICIT NONE

CONTAINS

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE compute_hg_adjoint ( Wd , id_Wd , QQ , u_opt , up_opt )

IMPLICIT NONE

TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(IN)  :: Wd
INTEGER                          , INTENT(IN)  :: id_Wd
TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(IN)  :: QQ
COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN)  :: u_opt
COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: up_opt

TYPE(CSR_MUMPS_Complex_Matrix) :: BB
COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: uRHS
INTEGER :: i1

ALLOCATE(uRHS  (velCmpnnts*np+np_L))
IF (ALLOCATED(up_opt)) DEALLOCATE(up_opt)
ALLOCATE(up_opt(SIZE(u_opt,1),SIZE(u_opt,2)))

DO i1 = 1 , SIZE(u_opt,2)
  CALL zZeroM ( QQ, np_L, BB)
  CALL zAtimx ( uRHS, BB%e, BB%j, BB%i, u_opt(:,i1) )
  uRHS = CONJG(uRHS) 
  CALL par_mumps_master(TRANSP_SOLUTION,id_Wd,Wd,0,uRHS)
  uRHS = CONJG(uRHS)
  up_opt(:,i1) = uRHS(1:velCmpnnts*np) 
ENDDO

END SUBROUTINE compute_hg_adjoint

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE compute_hgSens_baseFlow ( fflag , iom , Lns , QQ , id_QQ , &
                                     g_opt , f_opt , u_opt , up_opt )
! Some INPUTs may be unused depending on the case
! ...

IMPLICIT NONE

CHARACTER(*) , INTENT(IN) :: fflag
INTEGER      , INTENT(IN) :: iom
TYPE(CSR_MUMPS_Matrix)         , INTENT(IN) :: Lns
TYPE(CSR_MUMPS_Complex_Matrix) , INTENT(IN) :: QQ   
INTEGER                        , INTENT(IN) :: id_QQ
COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN) :: g_opt
COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: f_opt
COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: u_opt
COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: up_opt

TYPE(CSR_MUMPS_Complex_Matrix) :: Cu_c , Cu_v
COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: sens ! Sensitivity
                                                       ! to baseflow
                                                       ! modifications
COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: s1 , s2
COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: u_o
CHARACTER(LEN=24) :: format_string , omegaString
INTEGER :: n_hg
INTEGER :: i1

! number of HG
n_hg = p_in%hg_num_plot

! Allocatate Cu matrix ~ adjoint Oseen : C+(u*)
ALLOCATE(Cu_c%i      (SIZE(Lns%i))      ); Cu_c%i      =Lns%i
ALLOCATE(Cu_c%i_mumps(SIZE(Lns%i_mumps))); Cu_c%i_mumps=Lns%i_mumps
ALLOCATE(Cu_c%j      (SIZE(Lns%j))      ); Cu_c%j      =Lns%j
ALLOCATE(Cu_c%e      (SIZE(Lns%e))      )

ALLOCATE(u_o (velCmpnnts,np))
ALLOCATE(sens(velCmpnnts*np,n_hg)) ; sens = 0.0d0
ALLOCATE(s1(velCmpnnts*np),s2(velCmpnnts*np)) ; s1 = 0.0d0 ; s2 = 0.0d0

! Do loop on the HGs
DO i1 = 1 , n_hg

  CALL extract_z ( u_opt(:,i1) , u_o )

  ! No cumulation iter after iter :
  !  reinitialization is needed for FEM routines
  Cu_c%e = 0.0d0
  CALL qc_adjointOseen2_sp_M (mm,jj, &
             CMPLX(DBLE(u_o),-AIMAG(u_o),KIND=8),Cu_c)
  CALL zEssM ( Cu_c%e , Cu_c%j , Cu_c%i , velCmpnnts*np , &
               Cu_v%e , Cu_v%j , Cu_v%i , Cu_v%i_mumps  )
  
  IF     ( TRIM(fflag) == 'vol' ) THEN
    CALL zAtimx ( s1 , Cu_v%e , Cu_v%j , Cu_v%i , f_opt(:,i1) )
    s2 = - 2.0d0 * DBLE(g_opt(i1)) * DBLE(s1)
    CALL par_mumps_master (DIRECT_SOLUTION, id_QQ, QQ, 0, s2 )
    sens(:,i1) = DBLE(s2)
 
  ELSEIF ( TRIM(fflag) == 'in'  ) THEN
    CALL zAtimx ( s1 , Cu_v%e , Cu_v%j , Cu_v%i , up_opt(:,i1) )
    s2 = - 2.0d0 * DBLE(s1)
    CALL par_mumps_master (DIRECT_SOLUTION, id_QQ, QQ, 0, s2 )
    sens(:,i1) = DBLE(s2)

  ELSE
    WRITE(*,*) " ERROR in compute_hgSes_baseFlow   ! "
    WRITE(*,*) " fflag must be either 'vol' or 'in'  "
    STOP
  ENDIF

ENDDO

format_string = "(A2,I2.2)"
WRITE(omegaString,format_string) 'om' , iom 

CALL vtk_plot_eigenvectors ( rr , jj , DBLE(sens) ,             &
        TRIM(p_in%hg_output_directory) // TRIM(omegaString) //  &
                                                  'Usens.vtk' )

END SUBROUTINE compute_hgSens_baseFlow

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE compute_hgSens_steadyControl ( )

IMPLICIT NONE

WRITE(*,*) " compute_hgSens_steadyControl () TO BE WRITTEN YET "
STOP

END SUBROUTINE compute_hgSens_steadyControl

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE hg_compute_hg_sens
