MODULE qc_sp_M

!  COUPLED VELOCITY-PRESSURE SYSTEM OF EQUATIONS

!  THE USE OF THIS PROGRAM ASSUMES THAT THE LARGE SPARSE
!  MATRIX OF THE COUPLED SYSTEM IS STORED IN CSR FORMAT


!  ww(ni,l)      : Parabolic shape function
!  Dw_re(k,n,l)  : Derivative w.r.t x_k of ww on the reference simplex
!  d...          : Derivative in the physical space
!  pp_w(l)       : Weight of Gauss points of the Parabolic approximation
!  dwl(k,n)      : dw(n)/dxk * Jac_det(m)   ---> SUM(MNR(k,:,m)*Dw_re(:,n,l))
!
!  Don't forget to multiply by Jacobian determinant for mass terms
!  Don't forget to divide by Jacobian determinant for stiffness terms
!
!  For stiffness matrix only we have:
!  M^TM_j(k,h,m) = SUM(MNR(k,:,m)*MNR(h,:,m))/jac_det(m)
!
!   ===> dd_ij_ml = SUM_SUM(Dw_re(:,ni,l) * MTM_j(:,:,m) * Dw_re(:,nj,l))

  USE dynamic_structures

  USE sparse_matrix_profiles

 
CONTAINS

SUBROUTINE ComputeJacobianMatrix(np, mm, jj, jj_L, js_D, DESINGULARIZE, CC,  Re, UU)

   IMPLICIT NONE
   
   !-----------------------------------------------------------------------!
   INTEGER,                            INTENT(IN) :: np 
   INTEGER,            DIMENSION(:),   INTENT(IN) :: mm
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj_L
   TYPE(dyn_int_line), DIMENSION(:),   INTENT(IN) :: js_D 
   LOGICAL,                            INTENT(IN) :: DESINGULARIZE
   
   TYPE(CSR_MUMPS_Matrix),          INTENT(INOUT) :: CC
   
   REAL(KIND=8),                           INTENT(IN) :: Re
   REAL(KIND=8), OPTIONAL, DIMENSION(:,:), INTENT(IN) :: UU
   !-----------------------------------------------------------------------!

   REAL(KIND=8), PARAMETER :: zero = 0,  one = 1
   INTEGER :: Nx, p
   

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  [(U.G)_ + (_.G)U)]  + 1/Re * K_  +  G_ (weak)  
   !             CC  <---  - V._   
   !-------------ONLY THE CONSTANT CONTRIBUTION---------------------- 

   Nx = SIZE(CC%i) - 1

   CC%e = 0

   ! CALL qc_11_sp_M (mm, jj,     1d0/Re,  CC) ! + stiffness (ROT-DIV)

   CALL qc_11_sp_gg_M (mm, jj,  1d0/Re,  CC) ! + stifness (GRAD:GRAD)

   CALL qc_10_sp_M (mm, jj, jj_L, -one,  CC) ! + pressure gradient (ibp)

   CALL qc_01_sp_M (mm, jj, jj_L, -one,  CC) ! - velocity divergence
      
 
   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  [(U.G)_ + (_.G)U)]  + 1/Re * K_  +  G_ (ibp)  
   !-------------ADD THE ITERATION DEPENDENDT PART------------------- 


   IF (PRESENT(UU)) THEN
      CALL qc_oseen2_sp_M (mm, jj, UU,  CC) ! linearized terms
   END IF
   
   CALL Dirichlet_c_M (np, js_D,  CC)

   IF (DESINGULARIZE) THEN
      ! reduction of the row of the last equation to
      ! the diagonal element alone, set equal to 1
      DO p = CC%i(Nx), CC%i(Nx + 1) - 1
         CC%e(p) = 0
         IF (CC%j(p) == Nx) CC%e(p) = 1
      ENDDO
   ENDIF
   
END SUBROUTINE ComputeJacobianMatrix

!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!

SUBROUTINE ComputeJacobianMatrix_OLD(np, mm, jj, jj_L, js_D, DESINGULARIZE, CC,  Re_U)

   IMPLICIT NONE
   
   !-----------------------------------------------------------------------!
   INTEGER,                            INTENT(IN) :: np 
   INTEGER,            DIMENSION(:),   INTENT(IN) :: mm
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj_L
   TYPE(dyn_int_line), DIMENSION(:),   INTENT(IN) :: js_D 
   LOGICAL,                            INTENT(IN) :: DESINGULARIZE
   
   TYPE(CSR_MUMPS_Matrix),          INTENT(INOUT) :: CC
   
   REAL(KIND=8), OPTIONAL, DIMENSION(:,:), INTENT(IN) :: Re_U
   !-----------------------------------------------------------------------!

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: CCe_save
   REAL(KIND=8), PARAMETER :: zero = 0,  one = 1
   INTEGER :: Nx, p
   LOGICAL, SAVE :: initialized = .FALSE.
   

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  Re [(U.G)_ + (_.G)U)]  +  K_  +  G_ (weak)  
   !             CC  <---  - V._   
   !-------------ONLY THE CONSTANT CONTRIBUTION---------------------- 

   Nx = SIZE(CC%i) - 1

   IF (.NOT. initialized) THEN

      CC%e = zero  !  remind the incremental accumulation in CC
      !==========================================================

      initialized = .TRUE.

      CALL qc_11_sp_M (mm, jj,        one,  CC) ! + stiffness (ROT-DIV)

      CALL qc_10_sp_M (mm, jj, jj_L, -one,  CC) ! + pressure gradient (ibp)  

      CALL qc_01_sp_M (mm, jj, jj_L, -one,  CC) ! - velocity divergence
      
      ALLOCATE (CCe_save(SIZE(CC%e)))

      CCe_save = CC%e  ! STORE THE CONSTANT CONTRIBUTION (Stokes operator)
      
   END IF      

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  Re [(U.G)_ + (_.G)U)]  +  K_  +  G_ (ibp)  
   !-------------ADD THE ITERATION DEPENDENDT PART------------------- 

   CC%e = CCe_save

   IF (PRESENT(Re_U)) THEN 
      CALL qc_oseen2_sp_M (mm, jj, Re_U,  CC) ! linearized terms
   END IF
   
   CALL Dirichlet_c_M (np, js_D,  CC)

   IF (DESINGULARIZE) THEN    
      ! reduction of the row of the last equation to 
      ! the diagonal element alone, set equal to 1
      DO p = CC%i(Nx), CC%i(Nx + 1) - 1
         CC%e(p) = 0
         IF (CC%j(p) == Nx) CC%e(p) = 1
      ENDDO
   ENDIF
   
END SUBROUTINE ComputeJacobianMatrix_OLD

!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!

SUBROUTINE extract_Dirichlet_c (np, js_D, xx,  old_us_D)
!===========================================

   IMPLICIT NONE
   
   INTEGER,                           INTENT(IN)  :: np
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)  :: js_D
   REAL(KIND=8),        DIMENSION(:), INTENT(IN)  :: xx
   TYPE(dyn_real_line), DIMENSION(:)              :: old_us_D

   INTEGER :: k

   DO k = 1, 2

      old_us_D(k)%DRL = xx(js_D(k)%DIL + (k-1) * np)
   
   ENDDO
   
END SUBROUTINE extract_Dirichlet_c

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c_DIFF (np, js_D, us_D, old_us_D,  xx)
!===========================================

   IMPLICIT NONE
   
   INTEGER,                           INTENT(IN)  :: np
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)  :: js_D
   TYPE(dyn_real_line), DIMENSION(:), INTENT(IN)  :: us_D
   TYPE(dyn_real_line), DIMENSION(:), INTENT(IN)  :: old_us_D
   REAL(KIND=8),        DIMENSION(:)              :: xx

   INTEGER :: k

   DO k = 1, 2

      xx(js_D(k)%DIL + (k-1) * np) = us_D(k)%DRL - old_us_D(k)%DRL
   
   ENDDO
   
END SUBROUTINE Dirichlet_c_DIFF

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_homo_c_Cmplx (np, js_D,  xx)
!===========================================

   IMPLICIT NONE
   
   INTEGER,                           INTENT(IN)  :: np
   TYPE(dyn_int_line) , DIMENSION(:), INTENT(IN)  :: js_D
   COMPLEX(KIND=8)    , DIMENSION(:), INTENT(INOUT) :: xx

   INTEGER :: k

   DO k = 1, 2

      xx(js_D(k)%DIL + (k-1) * np) = 0.0d0 
  
   ENDDO
   
END SUBROUTINE Dirichlet_homo_c_Cmplx
!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c (np, js_D, us_D,  xx)
!===========================================

   IMPLICIT NONE
   
   INTEGER,                           INTENT(IN)  :: np
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)  :: js_D
   TYPE(dyn_real_line), DIMENSION(:), INTENT(IN)  :: us_D
   REAL(KIND=8),        DIMENSION(:), INTENT(OUT) :: xx

   INTEGER :: k

   DO k = 1, 2

      xx(js_D(k)%DIL + (k-1) * np) = us_D(k)%DRL
   
   ENDDO
   
END SUBROUTINE Dirichlet_c


!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c_M_MASS (np, js_D,  CC)
!=======================================

!  Modification of elements of selected rows of
!  matrix CC to impose Dirichlet boundary conditions 
!  at the nodes  js_D(1)%DIL  and  js_D(2)%DIL  

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,                           INTENT(IN)    :: np
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)    :: js_D
   TYPE(CSR_MUMPS_Matrix),            INTENT(INOUT) :: CC  

   INTEGER :: k, n, i_, p

  
   DO k = 1, 2
  
      DO n = 1, SIZE(js_D(k)%DIL)
    
         i_ = js_D(k)%DIL(n)  +  (k-1) * np
    
         DO p = CC%i(i_), CC%i(i_+1) - 1
         
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = 0d0
         
         ENDDO
    
      ENDDO
  
   ENDDO 
 
END SUBROUTINE Dirichlet_c_M_MASS

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c_M (np, js_D,  CC)
!=======================================

!  Modification of elements of selected rows of
!  matrix CC to impose Dirichlet boundary conditions 
!  at the nodes  js_D(1)%DIL  and  js_D(2)%DIL  

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,                           INTENT(IN)    :: np
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)    :: js_D
   TYPE(CSR_MUMPS_Matrix),            INTENT(INOUT) :: CC  

   INTEGER :: k, n, i_, p

  
   DO k = 1, 2
  
      DO n = 1, SIZE(js_D(k)%DIL)
    
         i_ = js_D(k)%DIL(n)  +  (k-1) * np
    
         DO p = CC%i(i_), CC%i(i_+1) - 1
         
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = 1
         
         ENDDO
    
      ENDDO
  
   ENDDO 
   
  
END SUBROUTINE Dirichlet_c_M

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_rc_M (np, js_D, diagE,  CC)
!=======================================

!  Modification of elements of selected rows AND columns
!  matrix CC to impose Dirichlet boundary conditions 
!  at the nodes  js_D(1)%DIL  and  js_D(2)%DIL  

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,                           INTENT(IN)    :: np
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)    :: js_D
   REAL(KIND=8),                      INTENT(IN)    :: diagE
   TYPE(CSR_MUMPS_Matrix),            INTENT(INOUT) :: CC 

   INTEGER :: k, n, i_, p

  
   DO k = 1, 2
  
      DO n = 1, SIZE(js_D(k)%DIL)
    
         i_ = js_D(k)%DIL(n)  +  (k-1) * np

         ! column
         WHERE ( CC%j == i_ )

            CC%e = 0

         ENDWHERE
    
         ! row
         DO p = CC%i(i_), CC%i(i_+1) - 1
         
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = diagE
         
         ENDDO
    
      ENDDO
  
   ENDDO 
   
  
END SUBROUTINE Dirichlet_rc_M

!------------------------------------------------------------------------------

SUBROUTINE qc_00_zero_sp_M (m0, jj, alpha,  CC) 
!==============================================


!  alpha << w, _ >>   ===>   CC  
!  
!  mass matrix but ONLY on Diagonal blocks 
!
!  The last (third) block which is zero (no cumulation)
!
!  ===>   CC   cumulative       

!  TWO-DIMENSIONAL VERSION ONLY

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: CC  
 
   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: a_w_w_p

   REAL(KIND=8) :: x
   INTEGER      :: np, mm, k, l, m, n, ni, nj, i, j, p, i_, j_


   IF (k_d /= 2) THEN
   
      WRITE (*,*) 'qc_00_zero_sp_M  is implemented only in 2D'
      WRITE (*,*) 'STOP.'
      STOP
   
   ENDIF


   np = MAXVAL(jj)  

   
   DO ni = 1, n_w

      DO nj = 1, n_w
      
         a_w_w_p(ni, nj, :) = alpha * ww(ni,:) * ww(nj,:) * pp_w
      
      ENDDO

   ENDDO


   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G
      
         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               x = a_w_w_p(ni, nj, l) * jac_det(m)

               ! diagonal block of the first block row
            
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
              
               ! diagonal block of the second blck row
               
               i_ = i + np;   j_ = j + np 
            
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_00_zero_sp_M

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE qc_11_sp_M (m0, jj, alpha,  CC) 
!=========================================

!  +  alpha [<< (Dxw), Dx_ >>  +  < (D.w), D._ >]   ===>   CC  
!                          
!  ===>   CC   cumulative       

!  TWO-DIMENSIONAL VERSION ONLY
!
!  +  alpha  times
!
!  dwx/dx . dvx/dx  +  dwx/dy . dvx/dy   |   dwx/dx . dvy/dy  -  dwx/dy . dvy/dx
!
! -dwy/dx . dvx/dy  +  dwy/dy . dvx/dx   |   dwy/dx . dvy/dx  +  dwy/dy . dvy/dy 
!
!
!  dwx/dx . d_x/dx  +  dwx/dy . d_x/dy   |   dwx/dx . d_y/dy  -  dwx/dy . d_y/dx
!
! -dwy/dx . d_x/dy  +  dwy/dy . d_x/dx   |   dwy/dx . d_y/dx  +  dwy/dy . d_y/dy 
!

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   
   REAL(KIND=8) :: alpha_pJ,  x, y
   INTEGER      :: np, mm, k, l, m, n, ni, nj, i, j, p, i_, j_

   IF (k_d /= 2) THEN
   
      WRITE (*,*) 'qc_11_sp_M  is implemented only in 2D'
      WRITE (*,*) 'STOP.'
      STOP
   
   ENDIF

   np = MAXVAL(jj)  
   
!   WRITE (*,*) 'np = ', np, 'Echo from qc_11_sp_M'
      

   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         alpha_pJ = alpha * pp_w(l) / jac_det(m)

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               ! FIRST BLOCK-ROW

               ! diagonal block 
               
               ! dwx/dx . dvx/dx  +  dwx/dy . dvx/dy
               ! dwx/dx . d_x/dx  +  dwx/dy . d_x/dy
               
               x = alpha_pJ * SUM(dwl(:,ni) * dwl(:,nj))

               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
               ! off-diagonal block of the first row  
                     
               j_ = j + np
              
               ! dwx/dx . dvy/dy  -  dwx/dy . dvy/dx
               ! dwx/dx . d_y/dy  -  dwx/dy . d_y/dx
             
               y = alpha_pJ * (dwl(1,ni) * dwl(2,nj) - dwl(2,ni) * dwl(1,nj))
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + y;  EXIT;  ENDIF
               ENDDO
               
               ! SECOND BLOCK-ROW
               
               i_ = i + np;
               
               ! off-diagonal block of the second row  
               
               ! -dwy/dx . dvx/dy  +  dwy/dy . dvx/dx
               ! -dwy/dx . d_x/dy  +  dwy/dy . d_x/dx
               
               ! x = alpha_pJ * (-dwl(1,ni) * dwl(2,nj) + dwl(2,ni) * dwl(1,nj))
               
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) - y;  EXIT;  ENDIF
               ENDDO                                     ! MINUS

               ! diagonal block 
              
               j_ = j + np 
              
               ! dwy/dx . dvy/dx  +  dwy/dy . dvy/dy 
               ! dwy/dx . d_y/dx  +  dwy/dy . d_y/dy 
             
               ! x = alpha_pJ * SUM(dwl(:,ni) * dwl(:,nj))
              
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_11_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_11_sp_gg_M (m0, jj, alpha,  CC) 
!===========================================
!
!  UNCOUPLED (DIAGONAL) STIFFNESS FOR THE VECTOR LAPLACIAN. 
!
!  IT CAN BE USED ONLY UNDER FULL DIRICHLET BOUNDARY CONDITIONS
!  FOR THE VECTOR UNKNOWN ALONG THE ENTIRE BOUNDARY

!  ACCUMULATES CONTRIBUTIONS ONLY TO DIAGONAL BLOCKS

!  +  alpha << (Dw), D_ >>   ===>   CC   /ALL VECTORS/
!                          
!  ===>   CC   cumulative       

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   
   REAL(KIND=8) :: x, alpha_pJ
   INTEGER      :: np, mm, k, l, m, n, ni, nj, i, j, p, i_, j_

   np = MAXVAL(jj)  
   
 !  WRITE (*,*) 'np = ', np, 'Echo from qc_11_sp_M'
      

   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         alpha_pJ = alpha * pp_w(l) / jac_det(m)

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               x = alpha_pJ * SUM(dwl(:,ni) * dwl(:,nj))

               ! first diagonal block 

               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
               ! second diagonal block 
               
               i_ = i + np;   j_ = j + np 
               
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_11_sp_gg_M

!------------------------------------------------------------------------------

SUBROUTINE qc_11_sp (m0, jj, gg, alpha,  v0)  !  BEWARE: NOT  _M
!====================================

!  + alpha [<< (Dxw), Dxg >>  +  < (D.w), D.g >]   ===>   v0     cumulative

!  TWO-DIMENSIONAL VERSION ONLY
!
!  +  alpha times
!
!  dwx/dx . dgx/dx  +  dwx/dy . dgx/dy   |   dwx/dx . dgy/dy  -  dwx/dy . dgy/dx
!
! -dwy/dx . dgx/dy  +  dwy/dy . dgx/dy   |   dwy/dx . dgy/dx  +  dwy/dy . dgy/dy 
!

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:,:)              :: v0
 
   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d, k_d) :: dgl
   REAL(KIND=8), DIMENSION(n_w)      :: dwdgl_k, dwdgl_k_OFF_DIAG
   INTEGER,      DIMENSION(n_w)      :: jjm
 
   INTEGER :: np, mm, k, k1, l, m, n
   
   IF (k_d /= 2) THEN
   
      WRITE (*,*) 'qc_11_sp  is implemented only in 2D'
      WRITE (*,*) 'STOP.'
      STOP
   
   ENDIF


   np = MAXVAL(jj)  
   
!   v0 = 0 

   DO mm = 1, SIZE(m0);  m = m0(mm)
      
      jjm = jj(:,m)

      ggm = gg(:,jjm)
  
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 1, k_d
            DO k1 = 1, k_d
               dgl(k,k1) = SUM(ggm(k,:) * dwl(k1,:)) 
            ENDDO
         ENDDO

         dgl = dgl * pp_w(l) / jac_det(m)  !  BEWARE 
     
         DO k = 1, k_d
          
            DO n = 1, n_w
            
               dwdgl_k(n) = SUM(dwl(:,n) * dgl(k,:))              
               
               SELECT CASE (k)
               
                  CASE(1);  dwdgl_k_OFF_DIAG(n) =  dwl(1,n) * dgl(2,2)  -  dwl(2,n) * dgl(2,1)
                  CASE(2);  dwdgl_k_OFF_DIAG(n) = -dwl(1,n) * dgl(1,2)  +  dwl(2,n) * dgl(1,1)
            
               END SELECT 
            
            ENDDO
         
            v0(k, jjm)  =  v0(k, jjm)  +  alpha * (dwdgl_k  +  dwdgl_k_OFF_DIAG)

         ENDDO
         
         ! dgl(1,1)  -->  dgx/dx  [times pp_w(l) / jac_det(m)]
         
         ! dgl(1,1)  -->  dgx/dx   |   dgl(1,2)  --->  dgx/dy
         ! dgl(2,1)  -->  dgy/dx   |   dgl(2,2)  --->  dgy/dy
         
         !  dw/dx . dgy/dy  -  dw/dy . dgy/dx
         ! -dw/dx . dgx/dy  +  dw/dy . dgx/dx        
         
         ! FIRST BLOCK-VECTOR
         !  dw/dx . dgx/dx  +  dw/dy . dgx/dy  |  dw/dx . dgy/dy  -  dw/dy . dgy/dx
 
         ! SECOND BLOCK-VECTOR   
         ! -dw/dx . dgx/dy  +  dw/dy . dgx/dx  |  dw/dx . dgy/dx  +  dw/dy . dgy/dy 
      
      ENDDO

   ENDDO

END SUBROUTINE qc_11_sp

!------------------------------------------------------------------------------

SUBROUTINE qc_11_sp_gg (m0, jj, gg, alpha,  v0)   !  BEWARE: NOT  _M 
!======================================   

!  UNCOUPLED (DIAGONAL) STIFFNESS FOR THE VECTOR LAPLACIAN. 
!
!  IT CAN BE USED ONLY UNDER FULL DIRICHLET BOUNDARY CONDITIONS
!  FOR THE VECTOR UNKNOWN ALONG THE ENTIRE BOUNDARY

!  alpha << (Dw), Dg >> 

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:,:)              :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl  
   REAL(KIND=8), DIMENSION(k_d, k_d) :: dgl
   REAL(KIND=8), DIMENSION(n_w)      :: dwdgl_k
   INTEGER,      DIMENSION(n_w)      :: jjm
 
   INTEGER :: mm, m, n, k, k1, l
      
!   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ggm = gg(:,jjm)
  
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 1, k_d
            DO k1 = 1, k_d
               dgl(k,k1) = SUM(ggm(k,:) * dwl(k1,:)) 
            ENDDO
         ENDDO
 
         DO k = 1, k_d
          
            DO n = 1, n_w
               dwdgl_k(n) = SUM(dwl(:,n) * dgl(k,:))              
            ENDDO
         
            v0(k, jjm)  =  v0(k, jjm)  + alpha * dwdgl_k * pp_w(l) / jac_det(m)

         ENDDO
         
      ENDDO

   ENDDO

END SUBROUTINE qc_11_sp_gg

!------------------------------------------------------------------------------

SUBROUTINE qc_advec_sp_M (m0, jj, gg,  CC) 
!=========================================

!  ACCUMULATES CONTRIBUTIONS ONLY TO DIAGONAL BLOCKS

!  +  << w, (g.D)_ >>   ===>   CC    /ALL VECTORS/
!                          
!  ===>   CC   cumulative      

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)      :: gl
   REAL(KIND=8), DIMENSION(n_w)      :: dgl_p
   INTEGER,      DIMENSION(n_w)      :: jjm

   INTEGER      :: np, mm, k, l, m, n, ni, nj, i, j, p, i_, j_
   REAL(KIND=8) :: x

   np = SIZE(gg, 2)

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO n = 1, n_w
            dgl_p(n) = SUM(gl * dwl(:,n)) * pp_w(l)
         ENDDO

         DO ni = 1, n_w;  i = jjm(ni)
                           
            DO nj = 1, n_w;  j = jjm(nj)

               x = ww(ni,l) * dgl_p(nj)
              
               ! first diagonal block  

               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

               ! second diagonal block 
               
               i_ = i + np;   j_ = j + np
                
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_advec_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_speci_sp_M (m0, jj, gg,  CC)
!=========================================

!  ACCUMULATES CONTRIBUTIONS TO ALL FOUR VELOCITY BLOCKS

!  +  << w, (_.D)g >>   ===>   CC    /ALL VECTORS/ 
!                          
!  ===>   CC   cumulative      

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl  
   REAL(KIND=8), DIMENSION(k_d, k_d) :: dgl_p
   INTEGER,      DIMENSION(n_w)      :: jjm

   INTEGER      :: np, mm, k, k1, l, m, n, ni, nj, i, j, p, i_, j_
   REAL(KIND=8) :: wij, x

   np = SIZE(gg, 2)

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 1, k_d
            DO k1 = 1, k_d
               dgl_p(k,k1) = SUM(dwl(k,:) * ggm(k1,:)) * pp_w(l)         
            ENDDO
         ENDDO
         
         DO ni = 1, n_w;  i = jjm(ni)
                           
            DO nj = 1, n_w;  j = jjm(nj)
               
               wij = ww(ni,l) * ww(nj,l)
               
               ! first diagonal block  
               
               ! x =  _ dgx/dx
               
               x = wij * dgl_p(1,1) 
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

               ! off diagonal block ux-uy of the first block-row  

               j_ = j + np

               ! x =  _ dgx/dy              

               x = wij * dgl_p(2,1)      
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
               ! off diagonal block uy-ux of the second first block-row 

               i_ = i + np 
               
               ! x =  _ dgy/dx   
             
               x = wij * dgl_p(1,2)     
             
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

               ! second diagonal block 

               j_ = j + np
            
               ! x =  _ dgy/dy                
             
               x = wij * dgl_p(2,2) 
             
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_speci_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_oseen2_sp_M (m0, jj, gg,  CC) 
!==========================================

!  ACCUMULATES CONTRIBUTIONS ONLY TO DIAGONAL BLOCKS

!  +  << w, (g.D)_ >>  +  << w, (_.D)g >>   ===>   CC    /ALL VECTORS/
!    
!  ACCUMULATES CONTRIBUTIONS TO ALL FOUR VELOCITY BLOCKS
!                      
!  ===>   CC   cumulative      

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl  
   REAL(KIND=8), DIMENSION(k_d, k_d) :: dgla_p
   REAL(KIND=8), DIMENSION(k_d)      :: gl
   REAL(KIND=8), DIMENSION(n_w)      :: dgl_p
   INTEGER,      DIMENSION(n_w)      :: jjm

   INTEGER      :: np, mm, k, k1, l, m, n, ni, nj, i, j, p, i_, j_
   REAL(KIND=8) :: wij, x, xa

   np = SIZE(gg, 2)

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
 
         DO n = 1, n_w
            dgl_p(n) = SUM(gl * dwl(:,n)) * pp_w(l)
         ENDDO

         DO k = 1, k_d
            DO k1 = 1, k_d
               dgla_p(k,k1) = SUM(dwl(k,:) * ggm(k1,:)) * pp_w(l)         
            ENDDO
         ENDDO
         
         DO ni = 1, n_w;  i = jjm(ni)
                           
            DO nj = 1, n_w;  j = jjm(nj)
               
               x = ww(ni,l) * dgl_p(nj) 
               
               wij = ww(ni,l) * ww(nj,l)
               
               ! first diagonal block  
               
               ! xa =  _ dgx/dx
                
               xa = wij * dgla_p(1,1) 
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x + xa;  EXIT;  ENDIF
               ENDDO

               ! off diagonal block ux-uy of the first block-row  

               j_ = j + np

               ! xa =  _ dgx/dy              

               xa = wij * dgla_p(2,1)      
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xa;  EXIT;  ENDIF
               ENDDO
               
               ! off diagonal block uy-ux of the second first block-row 

               i_ = i + np 
               
               ! xa =  _ dgy/dx   
             
               xa = wij * dgla_p(1,2)     
             
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + xa;  EXIT;  ENDIF
               ENDDO

               ! second diagonal block 

               j_ = j + np
            
               ! xa =  _ dgy/dy                
             
               xa = wij * dgla_p(2,2) 
             
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x + xa;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_oseen2_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_01P_sp_M (m0, jj, jj_L, alpha,  CC)
!================================================

!  +  alpha  << w, D_L >>   ===>   CC    
!                          
!  ===>   CC   cumulative      

   USE Gauss_points
   
   USE Gauss_points_L

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj, jj_L
   REAL(KIND=8),                 INTENT(IN)    :: alpha 
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(k_d, SIZE(jj_L,1), l_G) :: dw_L_re ! SIZE(jj_L,1) == n_w_L
   REAL(KIND=8), DIMENSION(k_d, SIZE(jj_L,1))      :: dw_Ll   ! SIZE(jj_L,1) == n_w_L
  
   INTEGER,      DIMENSION(n_w)          :: jjm
   INTEGER,      DIMENSION(SIZE(jj_L,1)) :: jjm_L             ! SIZE(jj_L,1) == n_w_L
  
   INTEGER      :: np, mm, m, l, n, k, ni, nj, i, j, p, i_, j_ 
   REAL(KIND=8) :: x
 
 
   np = MAXVAL(jj)  


   SELECT CASE (k_d)

      CASE (2)
               
         Dw_L_re(:,1,:) = Dw_re(:,1,:) + 0.5*(Dw_re(:,5,:) + Dw_re(:,6,:))
         Dw_L_re(:,2,:) = Dw_re(:,2,:) + 0.5*(Dw_re(:,6,:) + Dw_re(:,4,:))
         Dw_L_re(:,3,:) = Dw_re(:,3,:) + 0.5*(Dw_re(:,4,:) + Dw_re(:,5,:))
        
      !  w_L(1,:) = ww(1,:) + 0.5*(ww(5,:) + ww(6,:))
      !  w_L(2,:) = ww(2,:) + 0.5*(ww(6,:) + ww(4,:))
      !  w_L(3,:) = ww(3,:) + 0.5*(ww(4,:) + ww(5,:))

      CASE (3)
               
         Dw_L_re(:,1,:) = Dw_re(:,1,:) + 0.5*(Dw_re(:,n_w-2,:) + Dw_re(:,n_w-1,:) + Dw_re(:,n_w,:))
         Dw_L_re(:,2,:) = Dw_re(:,2,:) + 0.5*(Dw_re(:,6,:)     + Dw_re(:,n_w-3,:) + Dw_re(:,n_w,:))
         Dw_L_re(:,3,:) = Dw_re(:,3,:) + 0.5*(Dw_re(:,5,:)     + Dw_re(:,n_w-3,:) + Dw_re(:,n_w-1,:))
         Dw_L_re(:,4,:) = Dw_re(:,4,:) + 0.5*(Dw_re(:,5,:)     + Dw_re(:,6,:)     + Dw_re(:,n_w-2,:))

   END SELECT


   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm   = jj  (:,m)
      jjm_L = jj_L(:,m)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w_L !!!
               dw_Ll(k,n) = SUM(MNR(k,:,m) * Dw_L_re(:,n,l))              
            ENDDO
         ENDDO

         DO ni = 1, n_w;   i = jjm(ni);   i_ = i + np
              
            DO nj = 1, n_w_L;   j = jjm_L(nj);   j_ = j + 2*np

               ! first rectangular off-diagonal block of the last block-column  
                 
               x = alpha * ww(ni,l) * dw_Ll(1,nj) * pp_w(l)
                 
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
 
               ! second rectangular off-diagonal block of the last block-column  
 
               x = alpha * ww(ni,l) * dw_Ll(2,nj) * pp_w(l)
               
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_01P_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_10_sp_M (m0, jj, jj_L, alpha,  CC)
!===============================================

!  +  alpha < D.w, _L >   ===>   CC    
!                          
!  ===>   CC   cumulative      

   USE Gauss_points
   
   USE Gauss_points_L

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj, jj_L
   REAL(KIND=8),                 INTENT(IN)    :: alpha 
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(SIZE(jj_L,1), l_G) :: w_L

   REAL(KIND=8), DIMENSION(k_d, n_w)     :: dwl
   INTEGER,      DIMENSION(n_w)          :: jjm
   INTEGER,      DIMENSION(SIZE(jj_L,1)) :: jjm_L
  
   INTEGER      :: np, mm, m, l, n, k, ni, nj, i, j, p, i_, j_ 
   REAL(KIND=8) :: x
 
 
   np = MAXVAL(jj)  


   SELECT CASE (k_d)

      CASE (2)
               
         w_L(1,:) = ww(1,:) + 0.5*(ww(5,:) + ww(6,:))
         w_L(2,:) = ww(2,:) + 0.5*(ww(6,:) + ww(4,:))
         w_L(3,:) = ww(3,:) + 0.5*(ww(4,:) + ww(5,:))

      CASE (3)
               
         w_L(1,:) = ww(1,:) + 0.5*(ww(n_w-2,:) + ww(n_w-1,:) + ww(n_w,:))
         w_L(2,:) = ww(2,:) + 0.5*(ww(6,:)     + ww(n_w-3,:) + ww(n_w,:))
         w_L(3,:) = ww(3,:) + 0.5*(ww(5,:)     + ww(n_w-3,:) + ww(n_w-1,:))
         w_L(4,:) = ww(4,:) + 0.5*(ww(5,:)     + ww(6,:)     + ww(n_w-2,:))

   END SELECT


   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm   = jj  (:,m)
      jjm_L = jj_L(:,m)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO ni = 1, n_w;   i = jjm(ni);   i_ = i + np
              
            DO nj = 1, n_w_L;   j = jjm_L(nj);   j_ = j + 2*np

               ! first rectangular off-diagonal block of the last block-column  
                 
               x = alpha * dwl(1,ni) * w_L(nj,l) * pp_w(l)
                 
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
 
               ! second rectangular off-diagonal block of the last block-column  
 
               x = alpha * dwl(2,ni) * w_L(nj,l) * pp_w(l)
               
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_10_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_01_sp_M (m0, jj, jj_L, alpha,  CC)
!===============================================

!  +  alpha  < w_L, D._ >   ===>   CC   
!                          
!  ===>   CC   cumulative      

   USE Gauss_points
   
   USE Gauss_points_L

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj, jj_L
   REAL(KIND=8),                 INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(SIZE(jj_L,1), l_G) :: w_L

   REAL(KIND=8), DIMENSION(k_d, n_w)     :: dwl
   INTEGER,      DIMENSION(n_w)          :: jjm
   INTEGER,      DIMENSION(SIZE(jj_L,1)) :: jjm_L
  
   INTEGER      :: np, mm, m, l, n, k, ni, nj, i, j, p, i_, j_ 
   REAL(KIND=8) :: x  

   np = MAXVAL(jj)  
  

   SELECT CASE (k_d)

      CASE (2)
               
         w_L(1,:) = ww(1,:) + 0.5*(ww(5,:) + ww(6,:))
         w_L(2,:) = ww(2,:) + 0.5*(ww(6,:) + ww(4,:))
         w_L(3,:) = ww(3,:) + 0.5*(ww(4,:) + ww(5,:))

      CASE (3)
               
         w_L(1,:) = ww(1,:) + 0.5*(ww(n_w-2,:) + ww(n_w-1,:) + ww(n_w,:))
         w_L(2,:) = ww(2,:) + 0.5*(ww(6,:)     + ww(n_w-3,:) + ww(n_w,:))
         w_L(3,:) = ww(3,:) + 0.5*(ww(5,:)     + ww(n_w-3,:) + ww(n_w-1,:))
         w_L(4,:) = ww(4,:) + 0.5*(ww(5,:)     + ww(6,:)     + ww(n_w-2,:))

   END SELECT


   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm   = jj  (:,m)
      jjm_L = jj_L(:,m)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
 
         DO ni = 1, n_w_L;   i = jjm_L(ni);   i_ = i + 2*np

            DO nj = 1, n_w;   j = jjm(nj);   j_ = j + np
               
               ! first rectangular off-diagonal block of the bottom block-row  
              
               x = alpha * w_L(ni,l) * dwl(1,nj) * pp_w(l)
                 
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
  
               ! second rectangular off-diagonal block of the bottom block-row
                 
               x = alpha * w_L(ni,l) * dwl(2,nj) * pp_w(l)
                 
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_01_sp_M


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


SUBROUTINE qc_n0_sp_s (ms0, jjs, iis, fs,  v0)
!=============================================

!  << n.ws, fs >>_s   ===>   v0   incremental accumulation of boundary terms

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jjs, iis
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: fs
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: v0

   INTEGER :: mm, ms, ls, k
   REAL(KIND=8) :: fls

   DO mm = 1, SIZE(ms0);  ms = ms0(mm)
     
      DO ls = 1, l_Gs

         fls = SUM(fs(iis(:,ms)) * wws(:,ls)) * jac_dets(ms) * pp_ws(ls)

         DO k = 1, k_d     

            v0(k, jjs(:,ms)) = v0(k, jjs(:,ms)) + rnorms(k, ls,ms) * wws(:,ls) * fls

         ENDDO

      ENDDO
   
   ENDDO

       
END SUBROUTINE qc_n0_sp_s


!------------------------------------------------------------------------------


SUBROUTINE qc_t0_sp_s (ms0, jjs, iis, gzs,  v0)
!==============================================

!  << nxws, gzs >>_s   ===>   v0   incremental accumulation of boundary terms
  
   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jjs, iis
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: gzs
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8) :: gzls
   INTEGER :: mm, ms, ls
   
   IF (k_d /= 2) THEN
   
      WRITE (*,*) ' qc_t0_sp_s implemented only  in' 
      WRITE (*,*) ' the two-dimensional case: STOP.'
      STOP
       
   ENDIF
   
   DO mm = 1, SIZE(ms0);  ms = ms0(mm)

      DO ls = 1, l_Gs
         
         gzls = SUM(gzs(iis(:,ms)) * wws(:,ls)) * jac_dets(ms) * pp_ws(ls)
             
         v0(1, jjs(:,ms)) = v0(1, jjs(:,ms)) - rnorms(2, ls,ms) * wws(:,ls) * gzls
     
         v0(2, jjs(:,ms)) = v0(2, jjs(:,ms)) + rnorms(1, ls,ms) * wws(:,ls) * gzls
     
      ENDDO

   ENDDO


END SUBROUTINE qc_t0_sp_s

!==============================================================================

END MODULE qc_sp_M
