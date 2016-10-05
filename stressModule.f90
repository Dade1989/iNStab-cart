MODULE stressModule

 USE sparse_matrix_profiles
 USE global_variables
 USE prep_mesh_p1p2_sp
 USE Gauss_points
 USE Gauss_points_L
 USE qb_sp_m
 USE sparse_matrix_operations
 USE par_solve_mumps

IMPLICIT NONE

INTERFACE computeWallStress

  MODULE PROCEDURE computeWallStress_real  , &
                   computeWallStress_cmplx

END INTERFACE computeWallStress

INTERFACE computeForce

  MODULE PROCEDURE computeForce_real  , &
                   computeForce_cmplx

END INTERFACE computeForce


CONTAINS

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! computeWallStress_real ( xx , rr , side_Id , tauW , force ) 
! computeWallStress_cmplx( xx , rr , side_Id , tauW , force ) 
! computeForce_real      ( m0s , jjs , alpha , tauW , force )
! computeForce_cmplx     ( m0s , jjs , alpha , tauW , force )
! 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE computeWallStress_real ( xx , rr , side_Id , tauW , force)
  
  ! Real version
  
  IMPLICIT NONE
  
  REAL(KIND=8) , DIMENSION(:)   , INTENT(IN) :: xx
  REAL(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: rr
  INTEGER      , DIMENSION(:)   , INTENT(IN) :: side_Id
  
  REAL(KIND=8) , DIMENSION(:,:) , ALLOCATABLE:: tauW
  REAL(KIND=8) , DIMENSION(:)   , ALLOCATABLE:: force
  
  ! Test boundary connectivity ----------------------
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: bndEl
     INTEGER , DIMENSION(:)   , ALLOCATABLE :: iGlo
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: jjs5 , jj5Belt
     TYPE(CSR_MUMPS_Matrix) :: Mass_s , P_s , X1_s , X2_s ,X3_s , X4_s
     INTEGER :: i1 , k
  ! FEM problem -------------------------------------
     REAL(KIND=8) , DIMENSION(:) , ALLOCATABLE :: f1 , f2 , f3 , ff
  
  ! Check
  INTEGER :: i2
  INTEGER , DIMENSION(:) , ALLOCATABLE :: iV
  
  ! Test boundary connectivity and matrices----
  CALL volumeElementsOnSide(side_Id,jj,jjs,sides,bndEl,iGlo)
  
  ALLOCATE(jjs5(3,SIZE(bndEl,2)))
  jjs5 = bndEl(4:6,:)
  ALLOCATE(jj5Belt(6,SIZE(bndEl,2)))
  jj5Belt = jj(:,bndEl(3,:))
  
!  ! Check ----
!  OPEN(UNIT=20,FILE='./testRoutines/jj_jjs.txt')
!  WRITE(20,*) "        side         i1        mms(i1)          jjs(:,mms(i1))   "
!  DO i1 = 1 , SIZE(bndEl,2)
!    WRITE(20,*) bndEl(:,i1)
!  ENDDO
!  CLOSE(20)
!  ! ----------
  
  ! Start matrices (once for all)
  
  CALL start_matrix_P2bP2b ( jjs5 , Mass_s )
  
  CALL start_matrix_P2bP2belt ( jjs5 , jj5Belt , X1_s)
  CALL start_matrix_P2bP2belt ( jjs5 , jj5Belt , X2_s)
  CALL start_matrix_P2bP2belt ( jjs5 , jj5Belt , X3_s)
  CALL start_matrix_P2bP2belt ( jjs5 , jj5Belt , X4_s)
  
  CALL start_matrix_P2bP1b ( jjs5 , P_s )
  
  
  ! Matrix for the FEM problem:
  !  Mass matrix
  CALL qs_b_0b0b_sp_M  ( bndEl(2,:) , jjs ,     1.0d0 , Mass_s )
  
!  ! Check ----
!  OPEN(UNIT=20,FILE='./testRoutines/mass.txt')
!  WRITE(20,*) ' n_rows = ' , SIZE(Mass_s%i)-1
!  WRITE(20,*) ' nnz    = ' , SIZE(Mass_s%j)
!  WRITE(20,*) '       i         j        e     '
!  DO i1 = 1 , SIZE(Mass_s%j)
!    WRITE(20,*) Mass_s%i_mumps(i1) , Mass_s%j(i1) , Mass_s%e(i1)
!  ENDDO
!  CLOSE(20)
!  
!  ! ----------
  
  CALL par_mumps_master (INITIALIZATION, 10, Mass_s, 0)
  CALL par_mumps_master (SYMBO_FACTOR,   10, Mass_s, 0)
  CALL par_mumps_master (NUMER_FACTOR,   10, Mass_s, 0)
  
  !CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 2 , 2.0d0/Re , X1_s )
  !CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 2 , 1.0d0/Re , X1_s )
  
!  ! Check ----
!     X1_s%e = 0.0d0
!     CALL qs_b_0b1n_sp_M ( bndEl(3,:), jj5Belt, jjs5 , rr , 1 , 1 , 1.0d0 , X1_s )
!     X2_s%e = 0.0d0
!     CALL qs_b_0b1n_sp_M ( bndEl(3,:), jj5Belt, jjs5 , rr , 2 , 1 , 1.0d0 , X2_s )
!     X3_s%e = 0.0d0
!     CALL qs_b_0b1n_sp_M ( bndEl(3,:), jj5Belt, jjs5 , rr , 1 , 2 , 1.0d0 , X3_s )
!     X4_s%e = 0.0d0
!     CALL qs_b_0b1n_sp_M ( bndEl(3,:), jj5Belt, jjs5 , rr , 2 , 2 , 1.0d0 , X4_s )
!
!  Check ----
!  OPEN(UNIT=20,FILE='./testRoutines/aaadxnx.txt')
!  WRITE(20,*) ' n_rows = ' , SIZE(X1_s%i)-1
!  WRITE(20,*) ' nnz    = ' , SIZE(X1_s%j)
!  WRITE(20,*) '       i         j        e  d_dx ny  d_dy ny   u    v '
!  DO i1 = 1 , SIZE(X1_s%j)
!    WRITE(20,'( 2(I9) , 6(E17.5) )') X1_s%i_mumps(i1) , X1_s%j(i1) , & 
!                X1_s%e(i1),X2_s%e(i1) , &
!                X3_s%e(i1),X4_s%e(i1) , &
!                xx(X1_s%j(i1)) ,  xx(np+X1_s%j(i1))
!  ENDDO
!  CLOSE(20)
!  OPEN(UNIT=20,FILE='./testRoutines/dxnx40.txt')
!  WRITE(20,*) ' n_rows = ' , SIZE(X1_s%i)-1
!  WRITE(20,*) ' nnz    = ' , SIZE(X1_s%j)
!  WRITE(20,*) '       i         j        e  d_dx ny  d_dy ny   u    v '
!  ALLOCATE(iV(8))
!  iV = (/ 40 , 120 , 41 , 119 , 20 , 139 , 140 , 141 /)
!  DO i2 = 1 , 8
!  DO i1 = 1 , SIZE(X1_s%j)
!  IF ( X1_s%i_mumps(i1) .EQ. iV(i2) ) THEN
!    WRITE(20,'( 2(I9) , 6(E17.5) )') X1_s%i_mumps(i1) , X1_s%j(i1) , & 
!                X1_s%e(i1),X2_s%e(i1) , &
!                X3_s%e(i1),X4_s%e(i1) , &
!                xx(X1_s%j(i1)) ,  xx(np+X1_s%j(i1))
!  ENDIF
!  ENDDO
!  ENDDO
!  CLOSE(20)
!  ! ----------
  
  !STOP

  ! x component          ... ( ... , der , nCompnnt , ... ) 
  X1_s%e = 0.0d0
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 1 , 1 , 2.0d0/Re , X1_s )
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 2 , 1.0d0/Re , X1_s )
  X2_s%e = 0.0d0
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 1 , 2 , 1.0d0/Re , X2_s )
  
  P_s%e  = 0.0d0
  CALL qs_b_0b0bn_sp_M ( bndEl(2,:) , jjs , 1 , 1.0d0 , P_s )
  
  ! Build rhs_x
  ALLOCATE(f1(SIZE(Mass_s%i)-1))
  ALLOCATE(f2(SIZE(Mass_s%i)-1))
  ALLOCATE(f3(SIZE(Mass_s%i)-1))
  ALLOCATE(ff(SIZE(Mass_s%i)-1))
  
  
   CALL dAtimx ( f1 , X1_s%e , X1_s%j , X1_s%i , xx(1   :  np) )
   CALL dAtimx ( f2 , X2_s%e , X2_s%j , X2_s%i , xx(1+np:2*np) )
   CALL dAtimx ( f3 ,  P_s%e ,  P_s%j ,  P_s%i , xx(2*np+iGlo(1:MAXVAL(P_s%j)))  )
   ff = - f1 - f2 + f3
  ! ff = 1.0d0
  IF ( ALLOCATED(tauW) ) DEALLOCATE(tauW)
  ALLOCATE(tauW(2,SIZE(ff)))
  
  
  CALL par_mumps_master (DIRECT_SOLUTION, 10, Mass_s, 0, ff )
  tauW(1,:) = ff
  
  ! y component
  X1_s%e = 0.0d0
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 1 , 1.0d0/Re , X1_s )
  
  X2_s%e = 0.0d0
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 1 , 1 , 1.0d0/Re , X2_s )
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 2 , 2.0d0/Re , X2_s )
  
  P_s%e  = 0.0d0
  CALL qs_b_0b0bn_sp_M ( bndEl(2,:) , jjs , 2 , 1.0d0 , P_s )
  
   CALL dAtimx ( f1 , X1_s%e , X1_s%j , X1_s%i , xx(1   :  np) )
   CALL dAtimx ( f2 , X2_s%e , X2_s%j , X2_s%i , xx(1+np:2*np) )
   CALL dAtimx ( f3 ,  P_s%e ,  P_s%j ,  P_s%i , xx(2*np+iGlo(1:MAXVAL(P_s%j))) )
  
   ff = - f1 - f2 + f3
  !                 (/ (1.0d0 , k = 1 , MAXVAL(P_s%j) ) /)       )
  ! ff = 1.0d0
  
  CALL par_mumps_master (DIRECT_SOLUTION, 10, Mass_s, 0, ff) 
  tauW(2,:) = ff
  
  
  
  
!  ! Test -----
!  
!  WRITE(*,*)
!  WRITE(*,*) " Nodes 40 - 120 : s_x  , s_y "
!  WRITE(*,*) " i , xx , yy "
!  WRITE(*,*) iGlo(40)  , rr(1,iGlo(40))  , rr(2,iGlo(40))
!  WRITE(*,*) iGlo(120) , rr(1,iGlo(120)) , rr(2,iGlo(120))
!  WRITE(*,*) tauW(1,40 ) , tauW(2,40 )
!  WRITE(*,*) tauW(1,120) , tauW(2,120)
!  WRITE(*,*) " Nodes 41 - 119 : s_x  , s_y "
!  WRITE(*,*) iGlo(41)  , rr(1,iGlo(41)) , rr(2,iGlo(41))
!  WRITE(*,*) iGlo(119) , rr(1,iGlo(119)) , rr(2,iGlo(119))
!  WRITE(*,*) tauW(1,41 ) , tauW(2,41 )
!  WRITE(*,*) tauW(1,119) , tauW(2,119)
!  WRITE(*,*)
!  
!  
!  OPEN(UNIT=20,FILE='./testRoutines/stress.txt')
!  WRITE(20,*) "#   x   y     s_x        s_y       s_r       s_theta "
!  DO i1 = 1 , SIZE(ff)  
!    WRITE(20,*) rr(1,iGlo(i1)) , rr(2,iGlo(i1)) , tauW(1,i1) , tauW(2,i1) , &
!        (tauW(1,i1)*rr(1,iGlo(i1))+tauW(2,i1)*rr(2,iGlo(i1)))/ &
!        SQRT(rr(1,iGlo(i1))**2 + rr(2,iGlo(i1))**2 ) , &
!        (-tauW(1,i1)*rr(2,iGlo(i1))+tauW(2,i1)*rr(1,iGlo(i1)))/ &
!        SQRT(rr(1,iGlo(i1))**2 + rr(2,iGlo(i1))**2 )
!  ENDDO
!  CLOSE(20)
  
  IF (ALLOCATED(force)) DEALLOCATE(force) 
  ALLOCATE(force(2))
  
  CALL computeForce ( bndEl(2,:) , jjs5 , 1.0d0 , tauW(1,:) , force(1) )
  CALL computeForce ( bndEl(2,:) , jjs5 , 1.0d0 , tauW(2,:) , force(2) )
  
  CALL par_mumps_master (DEALLOCATION , 10 , Mass_s , 0 )
  
  ! ----------
  
  !   ! CHeck ----
  !   X1_s%e = 0.0d0
  !   CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 1 , 1 , 1.0d0 , X1_s )
  !   X2_s%e = 0.0d0
  !   CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 1 , 1.0d0 , X2_s )
  !   X3_s%e = 0.0d0
  !   CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 1 , 2 , 1.0d0 , X3_s )
  !   X4_s%e = 0.0d0
  !   CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 2 , 1.0d0 , X4_s )
  !   ! Check ----
  !   WRITE(*,*) " np , np_L , 2*np + np_L = " , np , np_L , 2*np+np_L
  !   WRITE(*,*) " SIZE(xx,1) =            = " , SIZE(xx,1)
  !   WRITE(*,*) " MINVAL(iGlo) , MAXVAL(iGlo) = " , MINVAL(iGlo) , MAXVAL(iGlo)
  !   WRITE(*,*) " MINVAL(iGlo(1:np_L)) , MAXVAL(iGlo(1:np_L)) = " ,  & 
  !                MINVAL(iGlo(1:MAXVAL(P_s%j))) , MAXVAL(iGlo(1:MAXVAL(P_s%j)))
  !   WRITE(*,*) " MINVAL(2*np+iGlo) , MAXVAL(2*np+iGlo) = " , &
  !                MINVAL(2*np+iGlo) , MAXVAL(2*np+iGlo)
  !   
  !   WRITE(*,*)  " MAXVAL(xx(1   :  np)) = " , MAXVAL(ABS(xx(1:np)))
  !   WRITE(*,*)  " MAXVAL(xx(1+np:2*np)) = " , MAXVAL(ABS(xx(1+np:2*np)))
  !   
  !   ! ----------
  
  
  
  END SUBROUTINE computeWallStress_real
  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  SUBROUTINE computeWallStress_cmplx ( xx , rr , side_Id , tauW , force)
  
  ! Complex version
  
  IMPLICIT NONE
  
  COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN) :: xx
  REAL(KIND=8)    , DIMENSION(:,:) , INTENT(IN) :: rr
  INTEGER         , DIMENSION(:)   , INTENT(IN) :: side_Id
  
  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE:: tauW
  COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE:: force
  
  ! Test boundary connectivity ----------------------
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: bndEl
     INTEGER , DIMENSION(:)   , ALLOCATABLE :: iGlo
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: jjs5 , jj5Belt
     TYPE(CSR_MUMPS_Matrix) :: Mass_s , P_s , X1_s , X2_s ,X3_s , X4_s
     TYPE(CSR_MUMPS_Complex_Matrix) :: Mass_c
     INTEGER :: i1 , k
  ! FEM problem -------------------------------------
     COMPLEX(KIND=8) , DIMENSION(:) , ALLOCATABLE :: f1 , f2 , f3 , ff
  
  ! Check
  INTEGER :: i2
  INTEGER , DIMENSION(:) , ALLOCATABLE :: iV
  
  ! Test boundary connectivity and matrices----
  CALL volumeElementsOnSide(side_Id,jj,jjs,sides,bndEl,iGlo)
  
  ALLOCATE(jjs5(3,SIZE(bndEl,2)))
  jjs5 = bndEl(4:6,:)
  ALLOCATE(jj5Belt(6,SIZE(bndEl,2)))
  jj5Belt = jj(:,bndEl(3,:))
  
!  ! Check ----
!  OPEN(UNIT=20,FILE='./testRoutines/jj_jjs.txt')
!  WRITE(20,*) "        side         i1        mms(i1)          jjs(:,mms(i1))   "
!  DO i1 = 1 , SIZE(bndEl,2)
!    WRITE(20,*) bndEl(:,i1)
!  ENDDO
!  CLOSE(20)
!  ! ----------
  
  ! Start matrices (once for all)
  
  CALL start_matrix_P2bP2b ( jjs5 , Mass_s )
  
  CALL start_matrix_P2bP2belt ( jjs5 , jj5Belt , X1_s)
  CALL start_matrix_P2bP2belt ( jjs5 , jj5Belt , X2_s)
  CALL start_matrix_P2bP2belt ( jjs5 , jj5Belt , X3_s)
  CALL start_matrix_P2bP2belt ( jjs5 , jj5Belt , X4_s)
  
  CALL start_matrix_P2bP1b ( jjs5 , P_s )
  
  
  ! Matrix for the FEM problem:
  !  Mass matrix
  CALL qs_b_0b0b_sp_M  ( bndEl(2,:) , jjs ,     1.0d0 , Mass_s )
  
  ALLOCATE(Mass_c%i      (SIZE(Mass_s%i      ))); Mass_c%i       = Mass_s%i
  ALLOCATE(Mass_c%i_mumps(SIZE(Mass_s%i_mumps))); Mass_c%i_mumps = Mass_s%i_mumps
  ALLOCATE(Mass_c%j      (SIZE(Mass_s%j      ))); Mass_c%j       = Mass_s%j
  ALLOCATE(Mass_c%e      (SIZE(Mass_s%e      )))
  Mass_c%e=CMPLX(Mass_s%e,0.0d0,KIND=8)

  CALL par_mumps_master (INITIALIZATION, 11, Mass_c, 0)
  CALL par_mumps_master (SYMBO_FACTOR,   11, Mass_c, 0)
  CALL par_mumps_master (NUMER_FACTOR,   11, Mass_c, 0)
  
! x component
! qs_b_0b1n_sp_M (m0s , jj , jjs , rr , deriv , nCmpnnt, alpha, AA)
  X1_s%e = 0.0d0
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 1 , 1, 2.0d0/Re , X1_s )
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 2, 1.0d0/Re , X1_s )
  X2_s%e = 0.0d0
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 1 , 2, 1.0d0/Re , X2_s )
  
  P_s%e  = 0.0d0
  CALL qs_b_0b0bn_sp_M ( bndEl(2,:) , jjs , 1 , 1.0d0 , P_s )
  
  ! Build rhs_x
  ALLOCATE(f1(SIZE(Mass_s%i)-1))
  ALLOCATE(f2(SIZE(Mass_s%i)-1))
  ALLOCATE(f3(SIZE(Mass_s%i)-1))
  ALLOCATE(ff(SIZE(Mass_s%i)-1))
  
  
  CALL zAtimx ( f1 , CMPLX(X1_s%e,0.0d0,KIND=8) , &
                           X1_s%j , X1_s%i, xx(1   :  np) )
  CALL zAtimx ( f2 , CMPLX(X2_s%e,0.0d0,KIND=8) , &
                           X2_s%j , X2_s%i, xx(1+np:2*np) )
  CALL zAtimx ( f3 , CMPLX( P_s%e,0.0d0,KIND=8) , &
                            P_s%j ,  P_s%i, xx(2*np+iGlo(1:MAXVAL(P_s%j)))  )
  ff = - f1 - f2 + f3
  ! ff = 1.0d0
  IF ( ALLOCATED(tauW) ) DEALLOCATE(tauW)
  ALLOCATE(tauW(2,SIZE(ff)))
  
  
  CALL par_mumps_master (DIRECT_SOLUTION, 11, Mass_c, 0, ff )
  tauW(1,:) = ff
  
! y component
  X1_s%e = 0.0d0
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 1, 1.0d0/Re , X1_s )
  
  X2_s%e = 0.0d0
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 1 , 1, 1.0d0/Re , X2_s )
  CALL qs_b_0b1n_sp_M ( bndEl(2,:), jj5Belt, jjs5 , rr , 2 , 2, 2.0d0/Re , X2_s )
  
  P_s%e  = 0.0d0
  CALL qs_b_0b0bn_sp_M ( bndEl(2,:) , jjs , 2 , 1.0d0 , P_s )
  
  CALL zAtimx ( f1 , CMPLX(X1_s%e,0.0d0,KIND=8) , &
                           X1_s%j , X1_s%i , xx(1   :  np) )
  CALL zAtimx ( f2 , CMPLX(X2_s%e,0.0d0,KIND=8) , &
                           X2_s%j , X2_s%i , xx(1+np:2*np) )
  CALL zAtimx ( f3 , CMPLX( P_s%e,0.0d0,KIND=8) , &
                            P_s%j ,  P_s%i , xx(2*np+iGlo(1:MAXVAL(P_s%j))) )
  
   ff = - f1 - f2 + f3
  !                 (/ (1.0d0 , k = 1 , MAXVAL(P_s%j) ) /)       )
  ! ff = 1.0d0
  
  CALL par_mumps_master (DIRECT_SOLUTION, 11, Mass_c, 0, ff) 
  tauW(2,:) = ff
  
  IF (ALLOCATED(force)) DEALLOCATE(force) 
  ALLOCATE(force(2))
  
  CALL computeForce ( bndEl(2,:) , jjs5 , CMPLX(1.0d0,0.0d0,KIND=8) , &
                                                       tauW(1,:) , force(1) )
  CALL computeForce ( bndEl(2,:) , jjs5 , CMPLX(1.0d0,0.0d0,KIND=8) , &
                                                       tauW(2,:) , force(2) )

  CALL par_mumps_master (DEALLOCATION , 11 , Mass_c , 0 )
  
  
  END SUBROUTINE computeWallStress_cmplx
  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  SUBROUTINE computeForce_real ( m0s , jjs , alpha , tauW , force )
  
  IMPLICIT NONE
  
  INTEGER      , DIMENSION(:)   , INTENT(IN) :: m0s
  INTEGER      , DIMENSION(:,:) , INTENT(IN) :: jjs
  REAL(KIND=8)                  , INTENT(IN) :: alpha
  REAL(KIND=8) , DIMENSION(:)   , INTENT(IN) :: tauW
  
  REAL(KIND=8) :: force
  
  
  CALL computeIntegral  ( m0s , jjs , alpha , tauW , force )
  
  
  END SUBROUTINE computeForce_real


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  SUBROUTINE computeForce_cmplx ( m0s , jjs , alpha , tauW , force )
  
  IMPLICIT NONE
  
  INTEGER         , DIMENSION(:)   , INTENT(IN) :: m0s
  INTEGER         , DIMENSION(:,:) , INTENT(IN) :: jjs
  COMPLEX(KIND=8)                  , INTENT(IN) :: alpha
  COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN) :: tauW
  
  COMPLEX(KIND=8) :: force
  
  
  CALL computeIntegral  ( m0s , jjs , alpha , tauW , force )
  
  
  END SUBROUTINE computeForce_cmplx

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

END MODULE stressModule
