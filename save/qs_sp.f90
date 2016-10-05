MODULE  qs_sp

CONTAINS

SUBROUTINE Dirichlet (js_D, us_D,  uu)
!=====================================

   IMPLICIT NONE
   INTEGER,      DIMENSION(:), INTENT(IN)  :: js_D
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: us_D
   REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: uu

   uu(js_D) = us_D

END SUBROUTINE Dirichlet

!------------------------------------------------------------------------------


SUBROUTINE qs_00_sp (m0, jj, ff,  u0)
!====================================

!  < w, f >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   REAL(KIND=8) :: fl
   INTEGER      :: mm, m, l

   u0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G

         fl = SUM(ff(jj(:,m)) * ww(:,l)) * jac_det(m) * pp_w(l)

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * fl

      ENDDO

   ENDDO

END SUBROUTINE qs_00_sp


!------------------------------------------------------------------------------


SUBROUTINE qs_01_sp (m0, jj, gg,  u0)
!====================================

!  < w, D.g >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl, ggm
   INTEGER,      DIMENSION(n_w)      :: jjm
  
   INTEGER      :: mm, m, l, n, k
   REAL(KIND=8) :: dgl_p

   u0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)
      
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         dgl_p = SUM(ggm * dwl) * pp_w(l)
         
         u0(jjm) = u0(jjm) + ww(:,l) * dgl_p
      
      ENDDO

   ENDDO


END SUBROUTINE qs_01_sp


!------------------------------------------------------------------------------


SUBROUTINE qs_11_sp (m0, jj, ff,  u0)
!====================================

!  < Dw, Df >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   REAL(KIND=8), DIMENSION(n_w)      :: ffm, dwdf_p  
   REAL(KIND=8), DIMENSION(k_d)      :: dfl
   INTEGER,      DIMENSION(n_w)      :: jjm
  
   INTEGER      :: mm, m, l, n, k

   u0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)
      
      ffm = ff(jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 1, k_d         
            dfl(k) = SUM(ffm * dwl(k,:))   
         ENDDO
 
         DO n = 1, n_w
            dwdf_p(n) = SUM(dwl(:,n) * dfl) * pp_w(l) / jac_det(m)
         ENDDO 
      
         u0(jjm) = u0(jjm) + dwdf_p 
      
      ENDDO

   ENDDO


END SUBROUTINE qs_11_sp


!------------------------------------------------------------------------------


SUBROUTINE qs_01_sp_c (m0, jj, gg,  u0)
!======================================

!  < w, D x g . k >>   ===>   u0       ( 2d only )

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl, ggm
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   INTEGER :: mm, m, l, k, n
   REAL(KIND=8) :: c_gl_p


   IF (k_d /= 2) THEN
   
       WRITE (*,*) 'Program qs_01_sp_c is valid only in two dimensions' 
   
       STOP 
       
   ENDIF


   u0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
    
      jjm = jj(:,m)
      
      ggm = gg(:, jjm)
    
      DO l = 1, l_G
    
         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
      
         c_gl_p = SUM(ggm(2,:) * dwl(1,:)   &
                    - ggm(1,:) * dwl(2,:)) * pp_w(l)

         u0(jjm) = u0(jjm) + ww(:,l) * c_gl_p

      ENDDO
 
   ENDDO


END SUBROUTINE qs_01_sp_c


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


SUBROUTINE qs_00_sp_s (ms0, jjs, iis, fs,  u0)
!=============================================

!  < ws, fs >_s   ===>   u0   incremental accumulation of boundary terms

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jjs, iis
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: fs
   REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: u0

   INTEGER :: mm, ms, ls
   REAL(KIND=8) :: fls

   DO mm = 1, SIZE(ms0);  ms = ms0(mm)
     
      DO ls = 1, l_Gs

         fls = SUM(fs(iis(:,ms)) * wws(:,ls)) * jac_dets(ms) * pp_ws(ls)

         u0(jjs(:,ms)) = u0(jjs(:,ms)) + wws(:,ls) * fls

      ENDDO
   
   ENDDO


END SUBROUTINE qs_00_sp_s


!------------------------------------------------------------------------------


SUBROUTINE qs_01_sp_s (ms0, jjs, iis, gs,  u0)
!=============================================

!  < ws, n.g_s >_s   ===>   u0   incremental accumulation of boundary terms
  
   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jjs, iis
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   REAL(KIND=8), DIMENSION(k_d) :: gls
   REAL(KIND=8) :: x
   INTEGER :: mm, ms, ls, k
   
   DO mm = 1, SIZE(ms0);  ms = ms0(mm)

      DO ls = 1, l_Gs
                                 
         DO k = 1, k_d
            gls(k) = SUM(gs(k, iis(:,ms)) * wws(:,ls))  
         ENDDO

         x = SUM(gls * rnorms(:,ls,ms)) * jac_dets(ms) * pp_ws(ls)
         
         u0(jjs(:,ms)) = u0(jjs(:,ms)) + wws(:,ls) * x 
                
      ENDDO

   ENDDO


END SUBROUTINE qs_01_sp_s



END MODULE qs_sp
