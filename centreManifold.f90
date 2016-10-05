MODULE centreManifold
  
 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Module to compute the centre manifold reduction for NS eqns:
!  the ROUTINES are (or at least should be) independent from
!  the dimension of the state and parameters.
! 
! It is known apriori that the critical space is composed of 
! two complex conjugate eigenvectors. One of them is stored in
! ./eigenOut/directEigenvectors, the other one is built as the
! c.c.
!  ( This module is first test on the Hopf bifurcation of the
!    flow around a circular cylinder , Re_c = 46.5...)
! BUT +++
!  Only very little changes (and only in reading the critical
!  space) should be required for higher-dimensionale critical
!  spaces.
! ++++++++++
!
! 2D CARTESIAN version
! cartesian and cylindrical version will differ only for the 
!  subroutine buildind the FEM matrices.
!
! 
! 2016-04-18 v1.0 : cylinder test ok
!  - cm_compute_cm_terms : resonant 01 term IF STATEMENT is missing
!
!     2016-04-19 : add stress computation and required routines
!      - add ROUTINEs to MODULE prep_mesh_p1p2_sp:
!         volumeElementsOnSide()
!      - add MODULE qb_sp_m: boundary FEM matrices
!
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! SUBROUTINE :
!  compute_centre_manifold()
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
  USE global_variables
  USE restart_io
  USE sparse_matrix_profiles
  USE par_solve_mumps

  USE cm_linear
  USE cm_forcing_input
  USE cm_index_connectivity
  USE cm_compute_cm_terms

  IMPLICIT NONE
 
  
 
  CONTAINS
  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  SUBROUTINE compute_centre_manifold()
  
  IMPLICIT NONE
  
  INTEGER :: nc
  
  TYPE(dyn_int_line),DIMENSION(velCmpnnts) :: js_D_eigen 
  TYPE(CSR_MUMPS_Complex_Matrix) :: Jacobian_cmplx, Mass_cmplx
  TYPE(CSR_MUMPS_Complex_Matrix) :: Ext_Jacobian_cmplx
  TYPE(CSR_MUMPS_Complex_Matrix) :: cm_Stiffness_cmplx
 
  INTEGER :: cmOrder 
  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: PhiC , PsiC
  COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: sigC
  
  INTEGER , DIMENSION(:)   , ALLOCATABLE :: a
  INTEGER , DIMENSION(:,:) , ALLOCATABLE :: bb
  INTEGER , PARAMETER :: npa = 1
  REAL(KIND=8) , DIMENSION(:) , ALLOCATABLE :: e0
  
  TYPE(powerExpansion) :: strIndTot
  INTEGER , DIMENSION(:,:) , ALLOCATABLE :: matIndTot
  
  INTEGER :: nm
  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: QQ , GG

  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: cmCF

  REAL(KIND=8) :: dummy
  INTEGER :: I1 , I2 , k , order
  
  ! End of Declarations
  
  
  ! Data needed to start Centre Manifold Reduction ++++++++++++++++++++
  !  Read centre manifold parameters : p_in % cm_max_order 
  cmOrder = p_in % cm_max_order
  
  !  Baseflow and parameters : baseflow = x0 , Re = Re
  CALL read_restart( x0, dummy, p_in%input_restart_file,  &
                     LEN(trim(p_in%input_restart_file)) )
  ALLOCATE(e0(npa))
  e0 = Re
  
  WRITE(*,*)  SIZE(x0,1) 
  WRITE(*,*) " e0 = " , e0
  !  Build or read the Mass and Jacobian Matrices
  CALL buildMassAndJacobian( x0 , js_D_eigen , Jacobian_Cmplx , Mass_Cmplx , & 
                     cm_Stiffness_Cmplx )

  !  Read (or compute) the critical subspace
  CALL readAndComputeCriticalSpace( PhiC , PsiC , sigC , nc )
  !  Eigenvectors normalization
  CALL eigenvectorNormalization ( Mass_cmplx%e , Mass_cmplx%j , &
                                  Mass_cmplx%i , PhiC , PsiC )

  !  Read nonlinear forcing terms
  CALL readForcingOrder ( a , bb ) 

 
  ! Centre Manifold Reduction +++++++++++++++++++++++++++++++++++++++++
  !  Connectivity
  CALL allTermsIndices  ( cmOrder , nc , npa , strIndTot )
  CALL fromMIstrToMImat ( strIndTot , matIndTot )
 
 
  !  Allocate structures
  nm = strIndTot % totTerms
  ALLOCATE(QQ(Nx,nm),GG(nc,nm))
  QQ = 0.0d0
  GG = 0.0d0
  
  WRITE(*,*)
  WRITE(*,*) " Centre Manifold Reduction : max order = " , cmOrder
  WRITE(*,*) " critical subspace         : n_crit    = " , nc
  WRITE(*,*) " parameter space           : n_par     = " , npa
  WRITE(*,*) " # terms of the expansion  : nm        = " , nm
  WRITE(*,*)

  ! MUMPS INITIALIZATIONS AND SYMBOLIC FACTORIZATIONs
  CALL par_mumps_master ( INITIALIZATION , 6 , Jacobian_cmplx, 0) 
  CALL par_mumps_master ( SYMBO_FACTOR   , 6 , Jacobian_cmplx, 0) 
  
  CALL start_extend_Jacobian  ( Jacobian_cmplx , 1 , Ext_Jacobian_Cmplx )

  CALL par_mumps_master ( INITIALIZATION , 7 , Ext_Jacobian_cmplx, 0) 
  CALL par_mumps_master ( SYMBO_FACTOR   , 7 , Ext_Jacobian_cmplx, 0) 
  
  ! order = 1
  WRITE(*,*)
  WRITE(*,*) " + Order = " , 1     , " ++++++++++++++++++++++++++++++++"
  !  (1,0)-terms from eigenpb
  DO i1 = 1 , nc ; QQ( :,i1) = PhiC(:,i1) ; ENDDO
  DO i1 = 1 , nc ; GG(i1,i1) = sigC(  i1) ; ENDDO
  
  !  (0,1)-terms
  CALL compute01terms ( x0 , e0 , matIndTot , js_D_eigen ,&
                      Jacobian_Cmplx , Ext_Jacobian_Cmplx , &
                      Mass_Cmplx , cm_Stiffness_Cmplx , &
                      PhiC , PsiC , sigC , &
                      a , bb , QQ , GG ) 
  
  DO order = 2 , cmOrder
    
    WRITE(*,*)
    WRITE(*,*) " + Order = " , order , " ++++++++++++++++++++++++++++++++"
    
    CALL computempterms ( x0 , e0 , order , matIndTot , js_D_eigen , &
                          Jacobian_Cmplx , Ext_Jacobian_Cmplx , &
                          Mass_Cmplx , cm_Stiffness_Cmplx , &
                          PhiC , PsiC , sigC , &
                          a , bb , QQ , GG )
    
  ENDDO


  OPEN(UNIT=22,FILE='./plots/GG.txt')
  WRITE(22,*) "       nc     ,       np      ,      nTerms  "
  WRITE(22,*) nc , npa ,SIZE(GG,2)
  WRITE(22,*) "        matIndTot                   GG"
  DO i1 = 1 , SIZE(GG,2)
    WRITE(22,*) matIndTot(i1,:) , GG(:,i1)
  ENDDO
  CLOSE(22)

  ! side_Id = 5 : cylinder
  CALL computeCFmpterms ( QQ , (/ 5 /) , cmCF )
   
  
  
END SUBROUTINE compute_centre_manifold

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE centreManifold
