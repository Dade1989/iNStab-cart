MODULE  qs_L_sp

CONTAINS

SUBROUTINE Dirichlet_L (js_D, us_D,  uu)
!=======================================

   IMPLICIT NONE

   INTEGER,      DIMENSION(:), INTENT(IN)  :: js_D
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: us_D
   REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: uu

   uu(js_D) = us_D 

END SUBROUTINE Dirichlet_L


!------------------------------------------------------------------------------


SUBROUTINE qs_01_hybrid_L_sp (m0, jj, jj_L, gg,  u0_L)
!=====================================================

!  < w_L, D.g >   ===>   u0_L

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj, jj_L
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:)                :: u0_L

   REAL(KIND=8), DIMENSION(SIZE(jj_L,1), l_G) :: w_L

   REAL(KIND=8), DIMENSION(k_d, n_w)     :: ggm, dwl
   REAL(KIND=8), DIMENSION(SIZE(jj_L,1)) :: uum_L
   INTEGER,      DIMENSION(n_w)          :: jjm
   INTEGER,      DIMENSION(SIZE(jj_L,1)) :: jjm_L
  
   INTEGER      :: mm, m, l, n, k, nwL, nw
   REAL(KIND=8) :: dgl_p

   nw  = SIZE(jj, 1)
   nwL = SIZE(jj_L, 1)


   SELECT CASE (k_d)


      CASE (2)
               
         w_L(1,:) = ww(1,:) + 0.5*(ww(5,:) + ww(6,:))
         w_L(2,:) = ww(2,:) + 0.5*(ww(6,:) + ww(4,:))
         w_L(3,:) = ww(3,:) + 0.5*(ww(4,:) + ww(5,:))

!        DO l = 1, l_G       
!           w_L(1,l) = ww(1,l) + 0.5*(ww(5,l) + ww(6,l))
!           w_L(2,l) = ww(2,l) + 0.5*(ww(6,l) + ww(4,l))
!           w_L(3,l) = ww(3,l) + 0.5*(ww(4,l) + ww(5,l))
!        END DO


      CASE (3)
               
         w_L(1,:) = ww(1,:) + 0.5*(ww(n_w-2,:) + ww(n_w-1,:) + ww(n_w,:))
         w_L(2,:) = ww(2,:) + 0.5*(ww(6,:)     + ww(n_w-3,:) + ww(n_w,:))
         w_L(3,:) = ww(3,:) + 0.5*(ww(5,:)     + ww(n_w-3,:) + ww(n_w-1,:))
         w_L(4,:) = ww(4,:) + 0.5*(ww(5,:)     + ww(6,:)     + ww(n_w-2,:))

!        DO l = 1, l_G     
!           w_L(1,l) = ww(1,l) + 0.5*(ww(n_w-2,l) + ww(n_w-1,l) + ww(n_w,l))
!           w_L(2,l) = ww(2,l) + 0.5*(ww(6,l)     + ww(n_w-3,l) + ww(n_w,l))
!           w_L(3,l) = ww(3,l) + 0.5*(ww(5,l)     + ww(n_w-3,l) + ww(n_w-1,l))
!           w_L(4,l) = ww(4,l) + 0.5*(ww(5,l)     + ww(6,l)     + ww(n_w-2,l))
!        END DO

  
   END SELECT


   !u0_L = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm   = jj  (:,m)
      jjm_L = jj_L(:,m)

      ggm = gg(:,jjm)

      uum_L = 0

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         dgl_p = SUM(ggm * dwl) * pp_w(l)

         DO n = 1, nwL
            uum_L(n) = uum_L(n) + w_L(n,l) * dgl_p
         ENDDO

      ENDDO

      u0_L(jjm_L) = u0_L(jjm_L) + uum_L

   ENDDO


END SUBROUTINE qs_01_hybrid_L_sp


!------------------------------------------------------------------------------


SUBROUTINE qs_00_L (m0, jj_L, ff_L,  u0_L)
!=========================================

!  < w_L, f_L >   ===>   u0_L

   USE Gauss_points_L

   IMPLICIT NONE
  
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_L
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff_L
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0_L
   
   REAL(KIND=8), DIMENSION(n_w_L) :: ffm_L
   INTEGER,      DIMENSION(n_w_L) :: jjm_L
   
   REAL(KIND=8) :: fl
   INTEGER :: mm, m, l
   
   u0_L = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
     
      jjm_L = jj_L(:,m)
     
      ffm_L = ff_L(jjm_L)
    
      DO l = 1, l_G_L

         fl = SUM(ffm_L * ww_L(:,l)) * rj_L(l,m)

         u0_L(jjm_L) = u0_L(jjm_L) + ww_L(:,l) * fl

      ENDDO
   
   ENDDO

END SUBROUTINE qs_00_L


END MODULE  qs_L_sp
