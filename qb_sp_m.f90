MODULE qb_sp_M

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Routines to compute FEM matrices built on boundary elements.
!
! Ex:
!   Matrix       |  weak form         |   discretized form
! ---------------+--------------------+--------------------------------
!  mass matrix   |  < w _  >_S        |  M_ij   =  < phi_i psi_j >_S
!  Gi            |  < w ni _ >_S      |  Gi_mn  =  < phi_m n_i psi_n >_S
!  Xij           |  < w nxj d_dxi >_S |  Xij_mn =  < phi_m nj dpsi_n/dxi >_S
! 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SUBROUTINEs :
!   countFromOne( A , B )
!   sort_pick   ( arr )
!- boundary matrices -----------
!   start_matrix_p2bp2b    ( jjswt , AA )
!   start_matrix_p2bp1b    ( jjswt , AA )
!- boundary - belt matrices ----
!   start_matrix_p2bp2belt ( jjswt , jj ,  AA )
!- boundary matrices -----------
!   qs_b_0b0b_sp_M    ( m0s ,      jjs ,                        alpha, AA)
!   qs_b_0b0bn_sp_M   ( m0s ,      jjs ,              nCmpnnt , alpha, AA)
!- boundary - belt matrices ----
!   qs_b_0b1n_sp_M    ( m0s , jj , jjs , rr , deriv , nCmpnnt , alpha, AA)
! ---  qs_b_0b0b_s_sp_M  ()   <----- boundary numbering ------------------
!                                    qc_00_WallOnly_sp_boundary_cmplx()
! ---                                in cylindrical code------------------
!   computeIntegral_real   ( m0s , jjs , alpha , v , AA ) 
!   computeIntegral_cmplx  ( m0s , jjs , alpha , v , AA ) 
!   localConnectivity ( jj , jjs , rr , jjN , MNRN , jac_detN )
! 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! The test functions phi are always defined on the boundary.
! The base functions psi could be the same as the test functions phi (when no
!  spatial derivative is required) or the base functions defined in the domain
!  (when spatial derivatives are required)
! Both the test and the base functions can be P1 or P2 or ...

  USE global_variables
  USE sparse_matrix_profiles
  
  IMPLICIT NONE
 
! Interfaces : computeIntegral is employed in stressModule to compute the
! integral actions on a wall or a whole body, e.g. lift or drag.
! If complex vector fields are integrated, the interface is required


  INTERFACE computeIntegral 
  
    MODULE PROCEDURE computeIntegral_real  , &
                     computeIntegral_cmplx
    
  END INTERFACE computeIntegral 



  
  CONTAINS

  
  SUBROUTINE countFromOne ( A , B )
     
     IMPLICIT NONE
     
     INTEGER , DIMENSION(:,:) , INTENT(IN)  :: A
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: B
     
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: ACOPY
     INTEGER :: I1 , I2 , I3
     INTEGER :: TOPV , MINV , IND , NUM
     
     IF ( ALLOCATED(B) ) THEN
        DEALLOCATE(B)
     END IF
     
     ALLOCATE(ACOPY(SIZE(A,1),SIZE(A,2)))
     ALLOCATE(B    (SIZE(A,1),SIZE(A,2)))
     B = 0
     ACOPY = A
     
     TOPV = MAXVAL(A) + 1
  !   MINV = MAXVAL(A) - 1
     IND  = 0
     NUM  = 0
     DO WHILE ( .NOT. (NUM == SIZE(ACOPY,1)*SIZE(ACOPY,2) ) )
        IND = IND + 1
        MINV = MINVAL(ACOPY)
        DO I2 = 1 , SIZE(ACOPY,1)
           DO I3 = 1 , SIZE(ACOPY,2)
              IF ( ACOPY(I2,I3) == MINV ) THEN
                 B(I2,I3)     = IND
                 ACOPY(I2,I3) = TOPV
                 NUM = NUM + 1
              ENDIF
           END DO
        END DO
     END DO
     
  
  END SUBROUTINE countFromOne  
  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE  sort_pick (arr)
  
  !  sorts an integer array  arr  into ascending numerical order, by stright
  !  insertion.   arr  is replaced on output by its sorted rearrangement.
  !  Not parallelizable.  Very inefficient for  SIZE(arr) > 20.
  
  !  Adapted from:  Numerical Recipes in Fortran 90,  p. 1167
  
     IMPLICIT NONE
  
     INTEGER, DIMENSION(:), INTENT(INOUT) :: arr
  
     INTEGER :: aj, i, j
  
     DO j = 2, SIZE(arr);   aj = arr(j)
  
        DO i = j - 1, 1, -1
           IF (arr(i) <= aj) EXIT
           arr(i+1) = arr(i)
        ENDDO
  
        arr(i+1) = aj
  
     ENDDO
  
  END SUBROUTINE  sort_pick

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE  start_matrix_p2bP2b ( JJSWT ,  AA )
  !+++++++++++++++++++++++++++++++++++++++++++++++
     
     IMPLICIT NONE
    
  !   INTEGER,                 INTENT(IN)  :: np
     INTEGER, DIMENSION(:,:)       , INTENT(IN)  :: JJSWT
     TYPE(CSR_MUMPS_Matrix), INTENT(OUT) :: AA
  !   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jjS
  !   INTEGER, DIMENSION(:),   INTENT(IN)  :: js
     
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: JJSW
     INTEGER :: NPp
     INTEGER, DIMENSION(:) , ALLOCATABLE :: nz_el_in_each_row
     INTEGER :: i, j, m, NA_nz, v, w, k, l, ii, c

     INTEGER :: i1 
     
     CALL countFromOne ( JJSWT , JJSW )
     
     npp = MAXVAL(JJSW)
     write(*,*) " # dofs on the boundary = " , npp
     
     ALLOCATE ( nz_el_in_each_row (npp) )
     ALLOCATE ( AA % i(npp+1) )
  
   
     ! preliminary determination of the boundary nodes that
     ! are mid-side nodes to have the right correction for
     ! all boundary nodes, including the mid-side ones. 
     
     nz_el_in_each_row = 1
! --- No preliminary determination of P1 node should be required:
! --- P1 nodes : 3 nz_el
! --- P2 nodes : 5 nz_el
! ---
!     DO m = 1, SIZE(jjSWT,2)
!        nz_el_in_each_row(jjsw(:,m))  =  nz_el_in_each_row(jjsw(:,m))  +  1
!     ENDDO
!     
!     WHERE (NZ_EL_IN_EACH_ROW .NE. 0)
!        nz_el_in_each_row = 1
!     ELSEWHERE
!        NZ_EL_IN_EACH_ROW = 0
!     END WHERE
! ---
! ---
     DO m = 1, SIZE(jjSWT,2)
  
        nz_el_in_each_row(jjSW(:,m))  =  nz_el_in_each_row(jjSW(:,m))  +  2
     
     ENDDO
     
  !   WRITE(*,*) NZ_EL_IN_EACH_ROW
     
     ! definition of the first (nz) element of each row
  
     AA % i(1) = 1
  
     DO i = 1, npp
  
        AA % i(i+1)  =  AA % i(i)  +  nz_el_in_each_row(i)
  
     ENDDO
  
     NA_nz  =  AA % i(npp+1)  -  1
!     write(*,*) "NA_nz = " , NA_nz
!     WRITE(*,*) "AA % i = "
!     WRITE(*,*) AA%i
  
     ALLOCATE (AA % j(NA_nz), AA % e(NA_nz))
     AA % j = 0
     AA % e = 0.0D0
  !+++ Edited by Pier   
     ALLOCATE ( AA % i_mumps (NA_nz) )
     
     c = 1
     
     DO i = 1, npp
     
        DO ii = AA%i(i), AA%i(i+1) -1
        
           AA%i_mumps(c) = i
           
           c = c + 1
        
        END DO
        
     END DO
  !+++
   
     !  column indices of all nonzero elements
  
  
     DO m = 1, SIZE(jjSWT,2)
  
        DO v = 1, SIZE(jjSWT,1)  ;   i = jjSW(v,m) !
  
           k  =  AA % i(i)           ! first element of the row
           l  =  AA % i(i+1)  -  1   ! last  element of the row
  
           DO w = 1 , SIZE(jjSW,1)  ;   J = jjSW(W,m)
  ! TEST ONLY ON P1 NODES ?????  ! SIZE(jj,1);   j = jj(w,m)
  !            write(*,*) 'jjsw(w,m) = ' , j
  !               write(*,*) 'AA % j(k:l) = ' , AA % j(k:l)
  !               write(*,*) 'jJSWT(W,M)  = ' , jJSWT(W,M)
              IF (ANY(AA % j(k:l) == jJSW(W,M)))   CYCLE
  
              ! forward shift and stack of the new element
  !            AA % j(l : k+1 : -1)  =  AA % j(l-1 : k : -1)
              AA % j(k+1 : l) = AA % j(k : l-1)
              AA % j(k)  =  j
  
           ENDDO
         
        ENDDO
  
     ENDDO
  
  
     ! sort in scending order of the column indices of each row
  
     DO i = 1, npp
  
        k  =  AA % i(i)           ! first element of the row
        l  =  AA % i(i+1)  -  1   ! last  element of the row
  
        CALL  sort_pick ( AA % j(k : l) )
  
     ENDDO
  
!     write (*,*) 'AA%j = ' , AA%j
     WRITE (*,*) '    CSR matrix structure for the P2S-P2s grid completed'

!   ! Test -----
!   OPEN(UNIT=20,FILE='./testRoutines/start_matrix_p2bp2b.txt')
!   WRITE(20,*) ' n_rows = ' , SIZE(AA%i)-1
!   WRITE(20,*) ' nnz    = ' , SIZE(AA%j)
!   WRITE(20,*) '       i         j       '
!   DO i1 = 1 , SIZE(AA%j)
!     WRITE(20,*) AA%i_mumps(i1) , AA%j(i1) 
!   ENDDO
!   CLOSE(20)
!   ! ----------

  
  
  END SUBROUTINE  start_matrix_p2bP2b

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE  start_matrix_p2bP1b ( JJSWT ,  AA )
  !+++++++++++++++++++++++++++++++++++++++++++++++
     
     IMPLICIT NONE
    
  !   INTEGER,                 INTENT(IN)  :: np
     INTEGER, DIMENSION(:,:)       , INTENT(IN)  :: JJSWT
     TYPE(CSR_MUMPS_Matrix), INTENT(OUT) :: AA
  !   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jjS
  !   INTEGER, DIMENSION(:),   INTENT(IN)  :: js
     
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: JJSW
     INTEGER :: NPp
     INTEGER, DIMENSION(:) , ALLOCATABLE :: nz_el_in_each_row
     INTEGER :: i, j, m, NA_nz, v, w, k, l, ii, c

     INTEGER :: i1 
     
     CALL countFromOne ( JJSWT , JJSW )
     
     npp = MAXVAL(JJSW)
     write(*,*) " # dofs on the boundary = " , npp
     
     ALLOCATE ( nz_el_in_each_row (npp) )
     ALLOCATE ( AA % i(npp+1) )
  
   
     ! preliminary determination of the boundary nodes that
     ! are mid-side nodes to have the right correction for
     ! all boundary nodes, including the mid-side ones. 
     
     nz_el_in_each_row = 1
! --- No preliminary determination of P1 node should be required:
! --- P1 nodes : 3 nz_el
! --- P2 nodes : 5 nz_el
! ---
!     DO m = 1, SIZE(jjSWT,2)
!        nz_el_in_each_row(jjsw(:,m))  =  nz_el_in_each_row(jjsw(:,m))  +  1
!     ENDDO
!     
!     WHERE (NZ_EL_IN_EACH_ROW .NE. 0)
!        nz_el_in_each_row = 1
!     ELSEWHERE
!        NZ_EL_IN_EACH_ROW = 0
!     END WHERE
! ---
! ---
     DO m = 1, SIZE(jjSWT,2)
  
        nz_el_in_each_row(jjSW(:,m))  =  nz_el_in_each_row(jjSW(:,m))  +  1
     
     ENDDO
     
  !   WRITE(*,*) NZ_EL_IN_EACH_ROW
     
     ! definition of the first (nz) element of each row
  
     AA % i(1) = 1
  
     DO i = 1, npp
  
        AA % i(i+1)  =  AA % i(i)  +  nz_el_in_each_row(i)
  
     ENDDO
  
     NA_nz  =  AA % i(npp+1)  -  1
!     write(*,*) "NA_nz = " , NA_nz
!     WRITE(*,*) "AA % i = "
!     WRITE(*,*) AA%i
  
     ALLOCATE (AA % j(NA_nz), AA % e(NA_nz))
     AA % j = 0
     AA % e = 0.0D0
  !+++ Edited by Pier   
     ALLOCATE ( AA % i_mumps (NA_nz) )
     
     c = 1
     
     DO i = 1, npp
     
        DO ii = AA%i(i), AA%i(i+1) -1
        
           AA%i_mumps(c) = i
           
           c = c + 1
        
        END DO
        
     END DO
  !+++
   
     !  column indices of all nonzero elements
  
  
     DO m = 1, SIZE(jjSWT,2)
  
        DO v = 1, SIZE(jjSWT,1)  ;   i = jjSW(v,m) !
  
           k  =  AA % i(i)           ! first element of the row
           l  =  AA % i(i+1)  -  1   ! last  element of the row
  
           DO w = 1 , 2             ;   J = jjSW(W,m)
  ! TEST ONLY ON P1 NODES ?????  ! SIZE(jj,1);   j = jj(w,m)
  !            write(*,*) 'jjsw(w,m) = ' , j
  !               write(*,*) 'AA % j(k:l) = ' , AA % j(k:l)
  !               write(*,*) 'jJSWT(W,M)  = ' , jJSWT(W,M)
              IF (ANY(AA % j(k:l) == jJSW(W,M)))   CYCLE
  
              ! forward shift and stack of the new element
  !            AA % j(l : k+1 : -1)  =  AA % j(l-1 : k : -1)
              AA % j(k+1 : l) = AA % j(k : l-1)
              AA % j(k)  =  j
  
           ENDDO
         
        ENDDO
  
     ENDDO
  
  
     ! sort in scending order of the column indices of each row
  
     DO i = 1, npp
  
        k  =  AA % i(i)           ! first element of the row
        l  =  AA % i(i+1)  -  1   ! last  element of the row
  
        CALL  sort_pick ( AA % j(k : l) )
  
     ENDDO
  
!     write (*,*) 'AA%j = ' , AA%j
     WRITE (*,*) '    CSR matrix structure for the P2S-P2s grid completed'
  
  
  END SUBROUTINE  start_matrix_p2bP1b

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE  start_matrix_p2bP2belt ( JJSWT , jj ,  AA )
  !+++++++++++++++++++++++++++++++++++++++++++++++
     
     IMPLICIT NONE
    
  !   INTEGER,                 INTENT(IN)  :: np
     INTEGER, DIMENSION(:,:) , INTENT(IN)  :: JJSWT
     INTEGER, DIMENSION(:,:) , INTENT(IN)  :: jj
     TYPE(CSR_MUMPS_Matrix)  , INTENT(OUT) :: AA
  !   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jjS
  !   INTEGER, DIMENSION(:),   INTENT(IN)  :: js
     
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: JJSW
     INTEGER :: NPp
     INTEGER, DIMENSION(:) , ALLOCATABLE :: nz_el_in_each_row
     INTEGER :: i, j, m, NA_nz, v, w, k, l, ii, c

     INTEGER :: i1 
     
     CALL countFromOne ( JJSWT , JJSW )

! Check ----
OPEN(UNIT=20,FILE='./testRoutines/jjsFromOne')
WRITE(20,*)
DO i1 = 1 , SIZE(jjswt,2)
  WRITE(20,*) i1
  WRITE(20,*) jjswt(:,i1) 
  WRITE(20,*) jjsw (:,i1) 
ENDDO
CLOSE(20)
! ----------

     
     npp = MAXVAL(JJSW)
     write(*,*) " # dofs on the boundary = " , npp
     
     ALLOCATE ( nz_el_in_each_row (npp) )
     ALLOCATE ( AA % i(npp+1) )
  
   
     ! preliminary determination of the boundary nodes that
     ! are mid-side nodes to have the right correction for
     ! all boundary nodes, including the mid-side ones. 
     
     nz_el_in_each_row = 1
! --- No preliminary determination of P1 node should be required:
! --- P1 nodes : 11 nz_el
! --- P2 nodes :  6 nz_el
! ---
!     DO m = 1, SIZE(jjSWT,2)
!        nz_el_in_each_row(jjsw(:,m))  =  nz_el_in_each_row(jjsw(:,m))  +  1
!     ENDDO
!     
!     WHERE (NZ_EL_IN_EACH_ROW .NE. 0)
!        nz_el_in_each_row = 1
!     ELSEWHERE
!        NZ_EL_IN_EACH_ROW = 0
!     END WHERE
! ---
! ---
     DO m = 1, SIZE(jjSWT,2)
  
        nz_el_in_each_row(jjSW(:,m))  =  nz_el_in_each_row(jjSW(:,m))  +  5
     
     ENDDO
     
  !   WRITE(*,*) NZ_EL_IN_EACH_ROW
     
     ! definition of the first (nz) element of each row
  
     AA % i(1) = 1
  
     DO i = 1, npp
  
        AA % i(i+1)  =  AA % i(i)  +  nz_el_in_each_row(i)
  
     ENDDO
  
     NA_nz  =  AA % i(npp+1)  -  1
!     write(*,*) "NA_nz = " , NA_nz
!     WRITE(*,*) "AA % i = "
!     WRITE(*,*) AA%i
  
     ALLOCATE (AA % j(NA_nz), AA % e(NA_nz))
     AA % j = 0
     AA % e = 0.0D0
  !+++ Edited by Pier   
     ALLOCATE ( AA % i_mumps (NA_nz) )
     
     c = 1
     
     DO i = 1, npp
     
        DO ii = AA%i(i), AA%i(i+1) -1
        
           AA%i_mumps(c) = i
           
           c = c + 1
        
        END DO
        
     END DO
  !+++
   
     !  column indices of all nonzero elements
  
  
     DO m = 1, SIZE(jjSWT,2)
  
        DO v = 1, SIZE(jjSWT,1)  ;   i = jjSW(v,m) !
  
           k  =  AA % i(i)           ! first element of the row
           l  =  AA % i(i+1)  -  1   ! last  element of the row
  
           DO w = 1 , SIZE(jj,1)    ;   J = jj(W,m)
  ! TEST ONLY ON P1 NODES ?????  ! SIZE(jj,1);   j = jj(w,m)
  !            write(*,*) 'jjsw(w,m) = ' , j
  !               write(*,*) 'AA % j(k:l) = ' , AA % j(k:l)
  !               write(*,*) 'jJSWT(W,M)  = ' , jJSWT(W,M)
              IF (ANY(AA % j(k:l) == jj(W,M)))   CYCLE
  
              ! forward shift and stack of the new element
  !            AA % j(l : k+1 : -1)  =  AA % j(l-1 : k : -1)
              AA % j(k+1 : l) = AA % j(k : l-1)
              AA % j(k)  =  j
  
           ENDDO
         
        ENDDO
  
     ENDDO
  
  
     ! sort in scending order of the column indices of each row
  
     DO i = 1, npp
  
        k  =  AA % i(i)           ! first element of the row
        l  =  AA % i(i+1)  -  1   ! last  element of the row
  
        CALL  sort_pick ( AA % j(k : l) )
  
     ENDDO
  
!     write (*,*) 'AA%j = ' , AA%j
     WRITE (*,*) '    CSR matrix structure for the P2S-P2s grid completed'
  
  
  END SUBROUTINE  start_matrix_p2bP2belt

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE qs_b_0b0b_sp_M (m0s, jjs, alpha, AA)
 
!!! It is assumed that the element in m0s are sorted in increasing order !!!

  ! mass matrix
  !  test functions : P2b
  !  base functions : P2b
  
  !  alpha < w, _ >_S    ===>   AA
  !                          
  !  ===>   AA   cumulative       
  
     USE Gauss_points
  
     IMPLICIT NONE
  
     INTEGER, DIMENSION(:),   INTENT(IN)    :: m0s
     INTEGER, DIMENSION(:,:), INTENT(IN)    :: jjs
     REAL(KIND=8),            INTENT(IN)    :: alpha
     TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: AA  

     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: jjs0
     INTEGER :: i0 , j0 

     REAL(KIND=8), DIMENSION(n_ws, n_ws, l_Gs) :: w_w_ps
  
     REAL(KIND=8) :: x, alpha_m
     INTEGER      :: mm, l, m, ni, nj, i, j, p
  

     CALL countFromOne( jjs(:,m0s) , jjs0 )
!  ! Check ----
!  WRITE(*,*) jjs(:,m0s)
!  WRITE(*,*) jjs0
!  ! ----------
     DO ni = 1, n_ws
  
        DO nj = 1, n_ws
        
           w_w_ps(ni, nj, :) = wws(ni,:) * wws(nj,:) * pp_ws
        
        ENDDO
  
     ENDDO

     DO mm = 1, SIZE(m0s);  m = m0s(mm)
!  ! Check ----
!  WRITE(*,*) m , mm , jac_dets(m)
!  ! ----------
        alpha_m = alpha * jac_dets(m)
  
        DO l = 1, l_Gs
  
           DO ni = 1, n_ws;  i = jjs(ni, m) ; i0 = jjs0(ni,mm)
  
              DO nj = 1, n_ws;  j = jjs(nj, m) ; j0 = jjs0(nj,mm)
  
                 x  =  alpha_m * w_w_ps(ni,nj,l)
  
                 DO p = AA%i(i0),  AA%i(i0+1) - 1
                   IF (AA%j(p) == j0) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
                 ENDDO
  
              ENDDO
  
           ENDDO
  
        ENDDO
  
     ENDDO
  
  
  END SUBROUTINE qs_b_0b0b_sp_M
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  SUBROUTINE qs_b_0b0bn_sp_M (m0s, jjs , nCmpnnt, alpha, AA)

!!! It is assumed that the element in m0s are sorted in increasing order !!!

  ! nCmpnnt : component of the normal unit vector in the integral
  
  !  
  !  test functions : P2b
  !  base functions : P1b
  
  !  alpha < w, ni _ >_S    ===>   AA
  !                          
  !  ===>   AA   cumulative       
  
     USE Gauss_points
     USE Gauss_points_L
  
     IMPLICIT NONE
  
     INTEGER, DIMENSION(:),   INTENT(IN)    :: m0s
     INTEGER, DIMENSION(:,:), INTENT(IN)    :: jjs
     INTEGER                , INTENT(IN)    :: nCmpnnt
     REAL(KIND=8),            INTENT(IN)    :: alpha
     TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: AA  
 
     REAL(KIND=8), DIMENSION(n_ws, n_ws_L, l_Gs) :: w_w_ps

     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: jjs0
     INTEGER :: i0 , j0 

     REAL(KIND=8) :: x, alpha_m
     INTEGER      :: mm, l, m, ni, nj, i, j, p
  
  !  To compute the P1 base fcns on the "P2" Gauss nodes
  !   REAL(KIND=8) , DIMENSION(SIZE(jjs_L,1), l_Gs) :: ws_L
     REAL(KIND=8) , DIMENSION(n_ws_L, l_Gs) :: ws_L
     REAL(KIND=8) , DIMENSION(l_Gs) :: xx
     REAL(KIND=8) :: f1 , f2 , y
 

      
  !  The value of the P1 base functions on the Gauss nodes for P2 are required
  ! P1 base functions : n_ws_L = 2 
     f1(y) = (1-y)/2.0d0 
     f2(y) = (1+y)/2.0d0 
  ! "P2" Gauss nodes : l_Gs = 3
     xx(1) = -SQRT(3.0d0/5.0d0)
     xx(2) =  0.0d0            
     xx(3) =  SQRT(3.0d0/5.0d0)
  
     DO ni = 1 , l_Gs
        ws_L(1,ni) = f1(xx(ni))
        ws_L(2,ni) = f2(xx(ni))
     ENDDO
  
  ! ----------------------------------------------------------------------------

     CALL countFromOne( jjs(:,m0s) , jjs0 )

     DO ni = 1, n_ws
  
        DO nj = 1, n_ws_L
        
           w_w_ps(ni, nj, :) = wws(ni,:) * ws_L(nj,:) * pp_ws
        
        ENDDO
  
     ENDDO
  
  
     DO mm = 1, SIZE(m0s);  m = m0s(mm)
  
        alpha_m = alpha * jac_dets(m)
  
        DO l = 1, l_Gs
  
           DO ni = 1, n_ws;  i = jjs(ni, m) ; i0 = jjs0(ni,mm)
  
              DO nj = 1, n_ws_L;  j = jjs(nj, m) ; j0 = jjs0(nj,mm)
  
                 x  =  alpha_m * rnorms(nCmpnnt,l,m) *  w_w_ps(ni,nj,l)
  
                 DO p = AA%i(i0),  AA%i(i0+1) - 1
                   IF (AA%j(p) == j0) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
                 ENDDO
  
              ENDDO
  
           ENDDO
  
        ENDDO
  
     ENDDO
  
  
  END SUBROUTINE qs_b_0b0bn_sp_M


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  SUBROUTINE qs_b_0b1n_sp_M (m0s , jj , jjs , rr , deriv , nCmpnnt, alpha, AA)

!!! It is assumed that the element in m0s are sorted in increasing order !!!

  ! m0s     : global numbering of the boundary elements on the desired side
  ! jj      : restricted to the BELT of domain elements with a side in common
  !           with the desired side
  ! jjs     : global numbering of the nodes on the desired side
  ! deriv   : coordinate of the partial derivative                : a
  ! nCmpnnt : component of the normal unit vector in the integral : b
  
  !  
  !  test functions : P2b
  !  base functions : P1b
  
  !  alpha < w, n_b d_/dx_a >_S    ===>   AA
  !                          
  !  ===>   AA   cumulative       
  
     USE Gauss_points
  
     IMPLICIT NONE
  
     INTEGER     , DIMENSION(:),    INTENT(IN) :: m0s
     INTEGER     , DIMENSION(:,:) , INTENT(IN) :: jj
     INTEGER     , DIMENSION(:,:) , INTENT(IN) :: jjs
     REAL(KIND=8), DIMENSION(:,:) , INTENT(IN) :: rr
     INTEGER                      , INTENT(IN) :: nCmpnnt
     INTEGER                      , INTENT(IN) :: deriv
     REAL(KIND=8),                  INTENT(IN) :: alpha
     TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: AA  
 
     REAL(KIND=8), DIMENSION(n_ws, n_ws, l_Gs) :: w_w_ps

     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: jjs0
     INTEGER :: i0 , j0 

     REAL(KIND=8) :: alpha_m , a
     INTEGER      :: mm, l, m, ni, nj, i, j, p , k

     INTEGER :: i1

! Intermediate connectivity
     INTEGER      , DIMENSION(:,:)   , ALLOCATABLE :: jjN
     REAL(KIND=8) , DIMENSION(:,:,:) , ALLOCATABLE :: MNRN
     REAL(KIND=8) , DIMENSION(:)     , ALLOCATABLE :: jac_detN
     REAL(KIND=8) , DIMENSION(:,:,:) , ALLOCATABLE :: dw_rg ! (k_d,n_w,n_Gs)
     REAL(KIND=8) , DIMENSION(:,:)   , ALLOCATABLE :: w_rg ! (k_d,n_w,n_Gs)

!  To compute the P1 base fcns on the "P2" Gauss nodes
     REAL(KIND=8) , DIMENSION(l_Gs) :: xx , yy , ee

     REAL(KIND=8) :: zero = 0,  half  = 0.5,  one  = 1,  &
                     two  = 2,  three = 3,    four = 4,  &
                     five = 5,   nine = 9

     REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,  &
                     df1x, df2x, df3x, df4x, df5x, df6x, &
                     df1y, df2y, df3y, df4y, df5y, df6y, &
                     x, y


! Base functions and derivatives for the domain elements : n_w
     f1(x, y) = (half - x - y) * (one - x - y) * two
     f2(x, y) = x * (x - half) * two
     f3(x, y) = y * (y - half) * two
     f4(x, y) = x * y * four
     f5(x, y) = y * (one - x - y) * four
     f6(x, y) = x * (one - x - y) * four

     df1x(x, y) = -three + four * (x + y)
     df2x(x, y) = (two*x - half) * two
     df3x(x, y) = zero
     df4x(x, y) =  y * four
     df5x(x, y) = -y * four
     df6x(x, y) = (one - two*x - y) * four

     df1y(x, y) = -three + four * (x + y)
     df2y(x, y) = zero
     df3y(x, y) = (two*y - half) * two
     df4y(x, y) =  x * four
     df5y(x, y) = (one - x - two*y) * four
     df6y(x, y) = -x * four


! "P2" Gauss nodes on the boundary : l_Gs = 3
     ee(1) = -SQRT(3.0d0/5.0d0)
     ee(2) =  0.0d0            
     ee(3) =  SQRT(3.0d0/5.0d0)

! The side of the physical element is brought to the side 12 of the
! reference element: 
! - 2d coordinates of the Gauss nodes on side 12
     DO i1 = 1 , 3
       xx(i1) = ( 1.0d0 + ee(i1) ) / two
       yy(i1) = 0.0d0
     ENDDO
! - Value of the domain base functions on the Gauss nodes on side 12
     ALLOCATE(dw_rg(2,6,3))   ! k_d,n_w,n_Gs
     ALLOCATE(w_rg(6,3))   ! k_d,n_w,n_Gs
     DO i1 = 1 , 3

      w_rg(1,i1) = f1(xx(i1),yy(i1))
      w_rg(2,i1) = f2(xx(i1),yy(i1))
      w_rg(3,i1) = f3(xx(i1),yy(i1))
      w_rg(4,i1) = f4(xx(i1),yy(i1))
      w_rg(5,i1) = f5(xx(i1),yy(i1))
      w_rg(6,i1) = f6(xx(i1),yy(i1))
     
      dw_rg(1,1,i1) = df1x(xx(i1),yy(i1)) ; dw_rg(2,1,i1) = df1y(xx(i1),yy(i1)) 
      dw_rg(1,2,i1) = df2x(xx(i1),yy(i1)) ; dw_rg(2,2,i1) = df2y(xx(i1),yy(i1)) 
      dw_rg(1,3,i1) = df3x(xx(i1),yy(i1)) ; dw_rg(2,3,i1) = df3y(xx(i1),yy(i1)) 
      dw_rg(1,4,i1) = df4x(xx(i1),yy(i1)) ; dw_rg(2,4,i1) = df4y(xx(i1),yy(i1)) 
      dw_rg(1,5,i1) = df5x(xx(i1),yy(i1)) ; dw_rg(2,5,i1) = df5y(xx(i1),yy(i1)) 
      dw_rg(1,6,i1) = df6x(xx(i1),yy(i1)) ; dw_rg(2,6,i1) = df6y(xx(i1),yy(i1)) 
     
     ENDDO


WRITE(*,*) MAXVAL(ABS(dw_rg))

! Intermediate connectivity 
     CALL localConnectivity( jj, jjs , rr , jjN , MNRN , jac_detN )

! Check ----
OPEN(UNIT=20,FILE='./testRoutines/MNRN.txt')
WRITE(20,*) 40 , jjN(:,40)
WRITE(20,*) 40 , MNRN(2,2,40) , -MNRN(1,2,40) , -MNRN(2,1,40) ,  MNRN(1,1,40) 
WRITE(20,*) 119 , jjN(:,119)
WRITE(20,*) 119,MNRN(2,2,119) ,-MNRN(1,2,119) ,-MNRN(2,1,119) ,  MNRN(1,1,119) 
WRITE(20,*) 39 , jjN(:,39)
WRITE(20,*) 39 , MNRN(2,2,39) ,-MNRN(1,2,39) ,-MNRN(2,1,39) ,  MNRN(1,1,39) 
WRITE(20,*) 120 , jjN(:,120)
WRITE(20,*) 120,MNRN(2,2,120) ,-MNRN(1,2,120) ,-MNRN(2,1,120) ,  MNRN(1,1,120) 


CLOSE(20)
! ----------

! Build the matrix : jjs already restricted to the desired boundary
     CALL countFromOne( jjs , jjs0 )


  DO mm = 1 , SIZE(m0s); m = m0s(mm)   ! sum over the elements
  
   DO l = 1 , l_Gs                      ! sum over the Gauss nodes
   
    DO k = 1 , 2                         ! sum over the partial derivatives
    
     DO ni = 1 , n_ws ; i = jjs(ni,mm) ; i0 = jjs0(ni,mm)   ! matrix rows
     
      DO nj = 1 , n_w  ; j = jj(nj,mm) ;  j0 = jjN(nj,mm)   ! matrix cols ??<--jj(nj,mm)
       
       a = alpha * pp_ws(l) * wws(ni,l) * & 
           dw_rg(k,nj,l) * &
           MNRN(deriv,k,mm)  / jac_detN(mm)  * &
           rnorms(nCmpnnt,l,m) * jac_dets(m)
 
       DO p = AA%i(i0) , AA%i(i0+1) - 1     ! cumulation
        IF (AA%j(p) == j0) THEN;  AA%e(p) = AA%e(p) + a;  EXIT;  ENDIF
       ENDDO
! Check ----
!       IF ( (mm .EQ. 37 .or. mm .EQ. 122) .AND. (i0 .EQ. 37 .or. i0.eq.123) ) THEN
!         WRITE(*,*) m
!         WRITE(*,*)  i0 , j0 , a
!WRITE(*,*) jac_dets(m) , jac_detN(mm)
!WRITE(*,*) rnorms(:,:,m) 
!WRITE(*,*) MNRN(:,:,mm) 
!       ENDIF
! ----------    
      ENDDO
     
     ENDDO
    
    ENDDO
   
   ENDDO
  
  ENDDO

 
 ! Test -----
 OPEN(UNIT=20,FILE='./testRoutines/boundary_domain.txt')
 WRITE(20,*) ' n_B_el    = ' , SIZE(jj,2)
 WRITE(20,*) ' jjs '
 WRITE(20,*) ' jj  '
 WRITE(20,*) ' jjN '
 DO i1 = 1 , SIZE(jj,2)
   WRITE(20,*)
   WRITE(20,*)  jjs(:,i1)
   WRITE(20,*)  jj (:,i1)
   WRITE(20,*)  jjN(:,i1)
   WRITE(20,*)  jac_det(m0s(i1)) , jac_detN(i1)
 ENDDO
 CLOSE(20)
 ! ----------
 
! ---------------------------------------------------------------------
  
  
  END SUBROUTINE qs_b_0b1n_sp_M

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE   qs_b_0b0b_s_sp_M (m0s, jj, jjs, alpha,  CC) 
!===============================================


!  alpha << w, _ >>|_{Gamma}   ===>   CC  
!  
!  mass matrix but ONLY on Diagonal blocks 
!

!  ===>   CC   cumulative       

!  TWO-DIMENSIONAL VERSION ONLY

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0s
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jjs
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_COMPLEX_Matrix),  INTENT(INOUT) :: CC  
   
   INTEGER, DIMENSION(:,:) , ALLOCATABLE :: JJSW
   INTEGER, DIMENSION(:,:) , ALLOCATABLE :: JJSWfromOne
   REAL(KIND=8), DIMENSION(n_ws, n_ws, l_Gs) :: a_w_w_p

   REAL(KIND=8) :: x
   INTEGER      :: npp, mm, l, m, ni, nj, i, j, p, i_, j_
   INTEGER      :: npRed
   
   INTEGER :: I1 , I2 , I3
   
   ALLOCATE(JJSW( SIZE(JJS,1), SIZE(M0S)) )
   JJSW = JJS(:,M0S)
   CALL countFromOne ( JJSW , JJSWfromOne )
   
   i2 = 0

   IF (k_d /= 2) THEN
   
      WRITE (*,*) 'qc_0y0_zero_sp_M  is implemented only in 2D'
      WRITE (*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   
   ENDIF


   npp    = MAXVAL(jj)  
   npRed  = MAXVAL(JJSWfromOne)
!   write(*,*) " npREd = " , npRed
   
   DO ni = 1, n_ws

      DO nj = 1, n_ws
      
         a_w_w_p(ni, nj, :) = alpha * wws(ni,:) * wws(nj,:) * pp_ws
      
      ENDDO

   ENDDO


   DO mm = 1, SIZE(m0s);  m = m0s(mm)            ! loop on the elements

      DO l = 1, l_Gs                             ! loop on the Gauss nodes
      
         DO ni = 1, n_ws;  i = jjsWFromOne(ni, mm)                                 ! loop on the element base function (1)

            DO nj = 1, n_ws;  j = jjsWFromOne(nj, mm)                               ! loop on the element base function (2)

               x = a_w_w_p(ni, nj, l) * jac_dets(m)    ! * yy_Gs(l,m)
!               x = alpha * wws(ni,l) * wws(nj,l) * pp_ws(l) * JACs(m)
               
               ! diagonal block of the first block row
            
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
              
               ! diagonal block of the second block row
               
               i_ = i + npRed;   j_ = j + npRed 
            
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
! 3d matrix   
!              ! diagonal block of the third block row
!              
!              i_ = i + 2*npRed;   j_ = j + 2*npRed
!           
!              DO p = CC%i(i_),  CC%i(i_+1) - 1
!                 IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
!              ENDDO
            
            ENDDO

         ENDDO

      ENDDO

   ENDDO


END SUBROUTINE  qs_b_0b0b_s_sp_M


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
  SUBROUTINE localConnectivity ( jj , jjs , rr , jjN , MNRN , jac_detN )
  
  ! jj
  ! jjs
  ! jjN      :  new local connectivity (required to make things simpler)
  ! MNRN     : 
  ! jac_detN :

  USE Gauss_points
 
  IMPLICIT NONE
  
  INTEGER      , DIMENSION(:,:)   , INTENT(IN)  :: jj , jjs
  REAL(KIND=8) , DIMENSION(:,:)   , INTENT(IN)  :: rr
  INTEGER      , DIMENSION(:,:)   , ALLOCATABLE :: jjN
  REAL(KIND=8) , DIMENSION(:,:,:) , ALLOCATABLE :: MNRN
  REAL(KIND=8) , DIMENSION(:)     , ALLOCATABLE :: jac_detN

  REAL(KIND=8) , DIMENSION(:,:)   , ALLOCATABLE :: drN

  INTEGER :: sz1 , sz2
  INTEGER :: ind
  INTEGER :: i1 , i2 , i3
  
  
  sz1 = SIZE(jj,1)
  sz2 = SIZE(jj,2)
  
  IF (ALLOCATED(jjN)) DEALLOCATE(jjN)
  ALLOCATE(jjN(sz1,sz2))
  jjN = 0 
  
  DO i1 = 1 , sz2
    
    DO i2 = 1 , 3 ! loop on the first 3 el of the jj
      
      DO i3 = 1 , 2 ! loop on the first 2 el of the jjs
        
        IF ( jj(i2,i1) .EQ. jjs(i3,i1) ) THEN
          jjN(i3  ,i1) = jj(i2  ,i1)
          jjN(i3+3,i1) = jj(i2+3,i1)
        ENDIF
        
        IF ( .NOT. ANY(jj(i2,i1).EQ.jjs(1:2,i1) ) ) THEN
          jjN(3  ,i1) = jj(i2  ,i1)
          jjN(3+3,i1) = jj(i2+3,i1)
        ENDIF
         
      ENDDO
      
    ENDDO
    
  ENDDO
  
! New MINOR matrix (MNRN) and determinants of the jacobian (det_jacN)
  ALLOCATE(MNRN(2,2,sz2)) ! k_d = 2 ; 2-dimensional pb
  ALLOCATE(jac_detN(sz2)) 
  ALLOCATE(drN(2,2)) 

  DO i1 = 1 , sz2

    drN(1,1) = rr(1,jjN(2,i1)) - rr(1,jjN(1,i1))
    drN(1,2) = rr(1,jjN(3,i1)) - rr(1,jjN(1,i1))
    drN(2,1) = rr(2,jjN(2,i1)) - rr(2,jjN(1,i1))
    drN(2,2) = rr(2,jjN(3,i1)) - rr(2,jjN(1,i1))

    MNRN(1,1,i1) =   drN(2,2)
    MNRN(1,2,i1) = - drN(2,1)
    MNRN(2,1,i1) = - drN(1,2)
    MNRN(2,2,i1) =   drN(1,1)

    jac_detN(i1) = drN(1,1) * drN(2,2) - drN(1,2) * drN(2,1)

  ENDDO 
  
  
  END SUBROUTINE localConnectivity

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE computeIntegral_real ( m0s , jjs , alpha , v , a )
  
  USE Gauss_points
  
  IMPLICIT NONE
  
  INTEGER      , DIMENSION(:)   , INTENT(IN) :: m0s
  INTEGER      , DIMENSION(:,:) , INTENT(IN) :: jjs
  REAL(KIND=8) ,                  INTENT(IN) :: alpha
  REAL(KIND=8) , DIMENSION(:)   , INTENT(IN) :: v
  
  REAL(KIND=8) :: a
  
  INTEGER      , DIMENSION(:,:) , ALLOCATABLE :: jjs0
  INTEGER :: mm , m , l , ni , p , i 
  
  a = 0.0d0
  
   CALL countFromOne( jjs , jjs0 )
  
   DO mm = 1 , SIZE(m0s); m = m0s(mm)   ! sum over the elements
    
     DO l = 1 , l_Gs                      ! sum over the Gauss nodes
      
       DO ni = 1 , n_ws ; i = jjs0(ni,mm)   ! matrix rows
         
         a = a + alpha * pp_ws(l) * wws(ni,l) * jac_dets(m) * v(i)
      
      ENDDO
     
     ENDDO
    
    ENDDO
  
  ! Check ----
  !       IF ( (mm .EQ. 37 .or. mm .EQ. 122) .AND. (i0 .EQ. 37 .or. i0.eq.123) ) THEN
  !         WRITE(*,*) m
  !         WRITE(*,*)  i0 , j0 , a
  !WRITE(*,*) jac_dets(m) , jac_detN(mm)
  !WRITE(*,*) rnorms(:,:,m) 
  !WRITE(*,*) MNRN(:,:,mm) 
  !       ENDIF
  ! ----------    
  
  
  END SUBROUTINE computeIntegral_real

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE computeIntegral_cmplx( m0s , jjs , alpha , v , a )
  
  USE Gauss_points
  
  IMPLICIT NONE
  
  INTEGER         , DIMENSION(:)   , INTENT(IN) :: m0s
  INTEGER         , DIMENSION(:,:) , INTENT(IN) :: jjs
  COMPLEX(KIND=8) ,                  INTENT(IN) :: alpha
  COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN) :: v
  
  COMPLEX(KIND=8) :: a
  
  INTEGER      , DIMENSION(:,:) , ALLOCATABLE :: jjs0
  INTEGER :: mm , m , l , ni , p , i 
  
  a = 0.0d0
  
   CALL countFromOne( jjs , jjs0 )
  
   DO mm = 1 , SIZE(m0s); m = m0s(mm)   ! sum over the elements
    
     DO l = 1 , l_Gs                      ! sum over the Gauss nodes
      
       DO ni = 1 , n_ws ; i = jjs0(ni,mm)   ! matrix rows
         
         a = a + alpha * pp_ws(l) * wws(ni,l) * jac_dets(m) * v(i)
      
      ENDDO
     
     ENDDO
    
    ENDDO
  
  END SUBROUTINE computeIntegral_cmplx

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE qb_sp_M
