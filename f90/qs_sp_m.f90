MODULE qs_sp_M

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

  USE sparse_matrix_profiles


  CONTAINS
 
 
SUBROUTINE Dirichlet_M (js_D, AA)
!================================

!  Mdification of elements of matrix AA to
!  impose Dirichlet boundary conditions at nodes js_D

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(IN)    :: js_D
   TYPE(CSR_MUMPS_Matrix),      INTENT(INOUT) :: AA  

   INTEGER :: n, i, p

   DO n = 1, SIZE(js_D);  i = js_D(n)
    
      DO p = AA%i(i), AA%i(i+1) - 1
        
         IF (AA%j(p) == i) THEN
            AA%e(p) = 1
         ELSE
            AA%e(p) = 0
         ENDIF
     
      ENDDO
   
   ENDDO

END SUBROUTINE Dirichlet_M


!------------------------------------------------------------------------------


SUBROUTINE qs_00_sp_M (m0, jj, alpha, AA)
!========================================

!  alpha < w, _ >    ===>   AA
!                          
!  ===>   AA   cumulative       

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: AA  

   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p

   REAL(KIND=8) :: x, alpha_m
   INTEGER      :: mm, l, m, ni, nj, i, j, p

   
   DO ni = 1, n_w

      DO nj = 1, n_w
      
         w_w_p(ni, nj, :) = ww(ni,:) * ww(nj,:) * pp_w
      
      ENDDO

   ENDDO


   DO mm = 1, SIZE(m0);  m = m0(mm)

      alpha_m = alpha * jac_det(m)

      DO l = 1, l_G

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               x  =  alpha_m * w_w_p(ni,nj,l)

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_00_sp_M


!------------------------------------------------------------------------------



SUBROUTINE qs_11_sp_M (m0, jj, alpha, AA)
!========================================

!  +  alpha << (Dw), D_ >>   ===>   AA
!                          
!  ===>   AA   cumulative       

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),        INTENT(INOUT) :: AA  

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   
   REAL(KIND=8) :: x, alpha_pj
   INTEGER      :: mm, k, l, m, n, ni, nj, i, j, p


   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         alpha_pj = alpha * pp_w(l) / jac_det(m)

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               x = alpha_pj * SUM(dwl(:,ni) * dwl(:,nj))

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_11_sp_M


!------------------------------------------------------------------------------




!------------------------------------------------------------------------------


SUBROUTINE qs_diff_mass_sp_M (m0, jj, beta, alpha, AA)
!=====================================================

!  +  beta << (Dw), D_ >>  +  alpha < w, _ >    ===>   AA
!                          
!  ===>   AA   cumulative       

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: beta, alpha
   TYPE(CSR_MUMPS_Matrix),        INTENT(INOUT) :: AA  

   REAL(KIND=8), DIMENSION(k_d, n_w)      :: dwl
   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p

   REAL(KIND=8) :: x, beta_p,  alpha_m,  dd_ij_ml
   INTEGER      :: mm, k, l, m, n, ni, nj, i, j, p

   
   DO ni = 1, n_w

      DO nj = 1, n_w
      
         w_w_p(ni, nj, :) = ww(ni,:) * ww(nj,:) * pp_w
      
      ENDDO

   ENDDO


   DO mm = 1, SIZE(m0);  m = m0(mm)

      alpha_m = alpha * jac_det(m)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO


         beta_p = beta * pp_w(l)

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               dd_ij_ml = SUM(dwl(:,ni) * dwl(:,nj))/jac_det(m)

               x  =  beta_p * dd_ij_ml  +  alpha_m * w_w_p(ni,nj,l)

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_diff_mass_sp_M


!------------------------------------------------------------------------------


SUBROUTINE qs_skew_adv_sp_M (m0, jj, gg, AA)
!===========================================

!  +  < w, (g.D)_ >  +  < w, (D.g) _ > / 2    ===>   AA
!                          
!  ===>   AA   cumulative      

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),             INTENT(INOUT) :: AA  

   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p
   REAL(KIND=8), DIMENSION(k_d, n_w)      :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)           :: gl
   REAL(KIND=8), DIMENSION(n_w)           :: gl_p
   INTEGER,      DIMENSION(n_w)           :: jjm

   INTEGER      :: mm, k, l, m, ni, nj, i, j, p, n
   REAL(KIND=8) :: dgl, x

   
   DO ni = 1, n_w
   
      DO nj = 1, n_w
         
         w_w_p(ni, nj, :) = ww(ni,:) * ww(nj,:) * pp_w / 2
               
      ENDDO
   
   ENDDO


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

         dgl = SUM(ggm * dwl)

         DO n = 1, n_w
            gl_p(n) = SUM(gl * dwl(:,n)) * pp_w(l)
         ENDDO

         DO ni = 1, n_w;  i = jjm(ni)

            DO nj = 1, n_w;  j = jjm(nj)

               x  =  ww(ni,l) * gl_p(nj)  +  dgl * w_w_p(ni,nj,l)

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_skew_adv_sp_M


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






END MODULE qs_sp_M




