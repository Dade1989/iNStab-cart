MODULE cm_linear

USE global_variables
USE sparse_matrix_profiles
USE Dirichlet_Neumann
USE qc_sp_M
USE cartesian_boundary_values
USE restart_io
USE eigensolve

IMPLICIT NONE

CONTAINS


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE start_extend_Jacobian ( jacobian , nPhi , ext_jacobian )

  IMPLICIT NONE
  
  TYPE(CSR_MUMPS_Complex_Matrix) , INTENT(IN) :: jacobian
  INTEGER                        , INTENT(IN) :: nPhi
  
  TYPE(CSR_MUMPS_Complex_Matrix)              :: ext_jacobian
  
  INTEGER :: nrows , i1 , k
  
  nrows = SIZE(jacobian%i)-1
  ALLOCATE( ext_jacobian%i      (SIZE(jacobian%i)      +nPhi        ) )
  ALLOCATE( ext_jacobian%i_mumps(SIZE(jacobian%i_mumps)+2*nPhi*nrows) )
  ALLOCATE( ext_jacobian%j      (SIZE(jacobian%j)      +2*nPhi*nrows) )
  ALLOCATE( ext_jacobian%e      (SIZE(jacobian%e)      +2*nPhi*nrows) )
  
  ! %i
  ext_Jacobian%i(1) = 1
  DO i1 = 1 , nrows
     ext_Jacobian%i(i1+1) = Jacobian%i(i1+1) + nPhi*i1
  ENDDO
  DO i1 = 1 , nPhi
     ext_Jacobian%i(nrows+1+i1) = ext_Jacobian%i(nrows+i1) + nRows
  ENDDO
 
 
  ! %i_mumps
  DO i1 = 1 , SIZE(ext_Jacobian%i) - 1
    ext_Jacobian%i_mumps(ext_Jacobian%i(i1):ext_Jacobian%i(i1+1)-1) = i1
  ENDDO
  
  ! %j
  DO i1 = 1 , nRows
    ext_Jacobian%j(ext_Jacobian%i(i1):ext_Jacobian%i(i1+1)-1-nPhi) =  &
        Jacobian%j(jacobian%i(i1):jacobian%i(i1+1)-1)
    ext_Jacobian%j(ext_Jacobian%i(i1+1)-nPhi:ext_Jacobian%i(i1+1)-1) = &
        (/ ( k , k = nRows+1 , nRows+nPhi) /)
  ENDDO
  DO i1 = 1 , nPhi
    ext_Jacobian%j(ext_Jacobian%i(nRows+i1):ext_Jacobian%i(nRows+i1+1)-1) = &
        (/ ( k , k = 1 , nRows ) /)
  ENDDO

  ! %e
  ext_Jacobian%e = 0.0d0

END SUBROUTINE start_extend_Jacobian


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE buildMassAndJacobian (x0,js_D_eigen, &
                    Jacobian_cmplx,Mass_cmplx,cm_Stiffness_cmplx)

  IMPLICIT NONE

  REAL(KIND=8) , DIMENSION(:) , INTENT(IN) :: x0
  TYPE(CSR_MUMPS_Complex_Matrix)           :: Jacobian_cmplx, Mass_cmplx
  TYPE(CSR_MUMPS_Complex_Matrix)           :: cm_Stiffness_cmplx
  LOGICAL, DIMENSION(velCmpnnts, number_of_sides) :: Dir_eigen
  TYPE(dyn_int_line), DIMENSION(velCmpnnts)       :: js_D_eigen
  LOGICAL                                         :: DESINGULARIZE_eigen
 
  COMPLEX(KIND=8) :: psum
  INTEGER :: k , i1 , i2
 
  ! Mass Matrix
  ! prepare boundary conditions
  !
  !   WRITE(*,*) '*check*'
  !   WRITE(*,*) '    number_of_sides = ', number_of_sides
  IF ( p_in%eigen_BC == 1 ) THEN
     ! homogeneous Dirichlet on every border
     WRITE(*,*) '    boundary conditions: zero velocity on every border'
     Dir_eigen = .TRUE.
     DESINGULARIZE_eigen = .TRUE.
  ELSEIF ( p_in%eigen_BC == 2 ) THEN
     ! same BCs as for the base flow but homogeneous
     ! Dirichlet where the base flow has nonhomogeneous Dirichlet
     WRITE(*,*) '    boundary conditions: same as base flow'
     Dir_eigen = Dir
     DESINGULARIZE_eigen = DESINGULARIZE
  ELSE
     WRITE(*,*) '*************************************'
     WRITE(*,*) '*** Wrong parameter:              ***'
     WRITE(*,*) '*** p_in % eigen_BC               ***'
     WRITE(*,*) '*** set to: ', p_in%eigen_BC
     WRITE(*,*) '*************************************'
     WRITE(*,*) 'STOP.'
     STOP
  ENDIF

  DO k = 1, velCmpnnts
    CALL Dirichlet_nodes_gen (jjs, sides, Dir_eigen(k,:), js_D_eigen(k)%DIL)
  ENDDO

  ! (2)
  ! update the Jacobian matrix with the solution received in input AND apply
  ! boundary conditions
  !
  Jacobian%e = 0
  CALL extract(x0, u0)

  WRITE(*,*) '*check*'
  WRITE(*,*) '    Re = ', Re
  CALL ComputeJacobianMatrix (np, mm, jj, jj_L, js_D_eigen,  &
                               DESINGULARIZE_eigen, Jacobian, Re, u0)
  CALL Dirichlet_rc_M (np, js_D_eigen, 1d0,  Jacobian)

  ALLOCATE( Jacobian_cmplx%i      (SIZE(Jacobian%i))       )
  ALLOCATE( Jacobian_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) )
  ALLOCATE( Jacobian_cmplx%j      (SIZE(Jacobian%j))       )
  ALLOCATE( Jacobian_cmplx%e      (SIZE(Jacobian%e))       )
  
  ALLOCATE( Mass_cmplx%i      (SIZE(Jacobian%i))       )
  ALLOCATE( Mass_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) )
  ALLOCATE( Mass_cmplx%j      (SIZE(Jacobian%j))       )
  ALLOCATE( Mass_cmplx%e      (SIZE(Jacobian%e))       )
  
  ALLOCATE( cm_Stiffness_cmplx%i      (SIZE(Jacobian%i))       )
  ALLOCATE( cm_Stiffness_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) )
  ALLOCATE( cm_Stiffness_cmplx%j      (SIZE(Jacobian%j))       )
  ALLOCATE( cm_Stiffness_cmplx%e      (SIZE(Jacobian%e))       )
  ! (4)
  ! convert the elements of the real Jacobian matrix (Jacobian) to Complex type
  ! copying them into the new CRS_MUMPS_Complex_Matrix Jacobian_cmplx AND
  ! change its sign as the eigenvalue problem we are solving is: 
  ! lambda*Mass*x = -Jacobian*x
  !
  Jacobian_cmplx%i       = Jacobian%i
  Jacobian_cmplx%i_mumps = Jacobian%i_mumps
  Jacobian_cmplx%j       = Jacobian%j
  Jacobian_cmplx%e       = CMPLX(-Jacobian%e, 0.0d0,KIND=8)
  
  ! (5)
  ! create the real Mass matrix with qc_0y0_zero_sp_M and
  ! convert the elements of the real Mass matrix (Mass) to Complex type
  ! copying them into the new CRS_MUMPS_Complex_Matrix Mass_cmplx
  !
  ! EXPLAINATION: we use the Jacobian matrix to save memory in case the Mass
  !               matrix hasn't been created yet
  !
  Jacobian%e = 0.0d0
  CALL qc_00_zero_sp_M (mm, jj, 1d0, Jacobian)
  
  ! impose boundary conditions on the Mass Matrix
  CALL Dirichlet_rc_M (np, js_D_eigen, 0d0,  Jacobian)
  
  Mass_cmplx%i       = Jacobian%i
  Mass_cmplx%i_mumps = Jacobian%i_mumps
  Mass_cmplx%j       = Jacobian%j
  Mass_cmplx%e       = CMPLX(Jacobian%e, 0d0, KIND=8)

  ! Stiffness Matrix
  Jacobian%e = 0.0d0
  CALL qc_11_sp_gg_M( mm , jj , 1d0 , Jacobian )

  ! impose boundary conditions on the Stiffness Matrix (NO b.c.)
!  CALL Dirichlet_c_M  (np, js_D_eigen,     Jacobian)
!  CALL Dirichlet_rc_M (np, js_D_eigen,0d0, Jacobian)

  cm_Stiffness_cmplx%i       = Jacobian%i
  cm_Stiffness_cmplx%i_mumps = Jacobian%i_mumps
  cm_Stiffness_cmplx%j       = Jacobian%j
  cm_Stiffness_cmplx%e       = CMPLX(Jacobian%e, 0d0, KIND=8)

!   ! Check ----
!   OPEN(UNIT=20,FILE="./plots/Ksum.txt")
!   WRITE(20,*) "sum of the elements on each row"
!   DO i1 = 1 , SIZE(cm_Stiffness_cmplx%i)-1
!   psum = 0.0d0
!   DO i2 = cm_Stiffness_cmplx%i(i1) , cm_Stiffness_cmplx%i(i1+1) -1 
!   psum = psum + cm_Stiffness_cmplx%e(i2)
!   ENDDO
!   WRITE(20,*) "row = " , i1 ," ; sum(row) = ", psum
!   ENDDO
!   CLOSE(20)
!   STOP
!   ! ----------
END SUBROUTINE buildMassAndJacobian

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE buildEyeMat ( L , Eye_Cmplx )

  IMPLICIT NONE

  TYPE(CSR_MUMPS_Complex_Matrix) , INTENT(IN)  :: L
  TYPE(CSR_MUMPS_Complex_Matrix)  :: Eye_cmplx
 
  INTEGER :: k , i1


  ALLOCATE( Eye_cmplx%i      (SIZE(L%i))      )
  ALLOCATE( Eye_cmplx%i_mumps(SIZE(L%i)-1)    )
  ALLOCATE( Eye_cmplx%j      (SIZE(L%i)-1)    )
  ALLOCATE( Eye_cmplx%e      (SIZE(L%i)-1)    )
  Eye_cmplx%i       = (/ ( k , k = 1,SIZE(L%i)  ) /)  
  Eye_cmplx%i_mumps = (/ ( k , k = 1,SIZE(L%i)-1) /)  
  Eye_cmplx%j       = Eye_cmplx%i_mumps 
  Eye_cmplx%e       = 1.0d0

! Check ---- 
! k = 0
! DO i1 = 1 , SIZE(Eye_Cmplx%e)
!   IF (Eye_Cmplx%e(i1) .EQ. 1.0d0 ) THEN ; k = k+1 ; END IF
! ENDDO
! WRITE(*,*) " rows of eye = " , SIZE(Eye_Cmplx%i) -1
! WRITE(*,*) " nnz el.     = " , k
! ----------

END SUBROUTINE buildEyeMat

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE readAndComputeCriticalSpace(PhiC,PsiC,sigC,nc)

  IMPLICIT NONE

  INTEGER :: nc  

  REAL(KIND=8) , DIMENSION(:) , POINTER :: rePhi1 
  REAL(KIND=8) , ALLOCATABLE , DIMENSION(:) , TARGET ::  imPhi1
  REAL(KIND=8) , DIMENSION(:) , ALLOCATABLE :: rePsi1 , imPsi1
  COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: sig1
  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: PhiC , PsiC
  COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: sigC

  ! Dimension of the critical subspace: if there is some zero 
  ! eigenvalue, this must be corrected
  ! nc = 2 * {non-zero critical eigenval} + {zero critical eigenval}
  nc = 2

  ALLOCATE(rePhi1(Nx),imPhi1(Nx))
  ALLOCATE(rePsi1(Nx),imPsi1(Nx))

  ! Read Critical Eigenspace (read)
  CALL read_eigenvector (Nx, 1, './eigenOut/Re46p51/directEigenvectors.dat', 41, &
                       rePhi1 , imPhi1 )
  CALL read_eigenvector (Nx, 1, './eigenOut/Re46p51/adjointEigenvectors.dat', 42,&
                       rePsi1 , imPsi1 )
  CALL read_eigenvalues ('./eigenOut/Re46p51/directEigenvalues.dat',1,sig1)
  

  ! PhiC = [phi , CONJG(phi)]
  ALLOCATE(PhiC(SIZE(rePhi1),2),PsiC(SIZE(rePsi1),2))
  ALLOCATE(sigC(2))
  PhiC(:,1) = CMPLX(rePhi1,imPhi1) ; PhiC(:,2) = CMPLX(rePhi1,-imPhi1) 
  PsiC(:,1) = CMPLX(rePsi1,imPsi1) ; PsiC(:,2) = CMPLX(rePsi1,-imPsi1)
  sigC(1) = sig1(1) ; sigC(2) = CONJG(sig1(1))
  
  WRITE(*,*)
  WRITE(*,*) " Critical eigenvalues = "
  WRITE(*,*) sigC(1)
  WRITE(*,*) sigC(2)
  WRITE(*,*)

END SUBROUTINE readAndComputeCriticalSpace

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE eigenvectorNormalization ( b , bj , bi , Phi , Psi )

  IMPLICIT NONE
  
  COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN) :: b
  INTEGER         , DIMENSION(:)   , INTENT(IN) :: bi , bj
  COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(INOUT) :: PHI , PSI
  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: BPhi
  
  INTEGER :: nval
  
  INTEGER :: i1
  
  nval = SIZE(PHI,2)
  
  ALLOCATE(BPHI(SIZE(Phi,1),SIZE(Phi,2)) )
  BPhi = 0.0d0

  ! Normalization of the right eigenvector |Phi_j| = 1
  DO i1 = 1 , nval
    CALL zAtimx ( BPhi(:,I1) , b , bj , bi , CONJG(PHI(:,I1)) )
    PHI(:,I1) = PHI(:,I1) / SQRT( SUM( PHI(:,I1) * BPhi(:,i1) ) )
  ENDDO
  
  WRITE(*,*) SUM(PHI(:,1) * BPHI(:,1)) , SUM(PHI(:,2) * BPHI(:,1))
  WRITE(*,*) SUM(PHI(:,1) * BPHI(:,2)) , SUM(PHI(:,2) * BPHI(:,2))

  ! Normalization to unitary mass
  DO I1 = 1 , NVAL
    CALL zAtimx ( BPhi(:,I1) , b , bj , bi , CONJG(PHI(:,I1)) )
    PSI(:,I1) = PSI(:,I1) / ( SUM(PSI(:,I1) * BPhi(:,I1) ) )
  ENDDO 

  ! Check orthogonality
  
  WRITE(*,*) SUM(PSI(:,1) * BPHI(:,1)) , SUM(PSI(:,2) * BPHI(:,1))
  WRITE(*,*) SUM(PSI(:,1) * BPHI(:,2)) , SUM(PSI(:,2) * BPHI(:,2))

END SUBROUTINE eigenvectorNormalization

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE cm_linear
