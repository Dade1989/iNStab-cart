MODULE qs_L_M

   
   USE sparse_matrix_profiles
   
   
   CONTAINS


SUBROUTINE Dirichlet_L_M (js_D, AA)
!==================================

!  Enforces Dirichlet boundary conditions
   
   IMPLICIT NONE

   INTEGER, DIMENSION(:),  INTENT(IN)    :: js_D
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA  

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
  
END SUBROUTINE Dirichlet_L_M


!------------------------------------------------------------------------------


SUBROUTINE qs_00_L_M (m0, jj, alpha, AA)
!=======================================

!  alpha < w, _ >   ===>   AA    incremental accumulation

   USE Gauss_points
   USE Gauss_points_L

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),                 INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: AA
   
   INTEGER :: mm, m, l, ni, nj, i, j, p
   REAL(KIND=8) :: al, x

   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G_L

         al = alpha * rj_L(l,m)

         DO ni = 1, n_w_L;  i = jj(ni, m)

            DO nj = 1, n_w_L;  j = jj(nj, m)

!              IF (j >= i) THEN
             
               x = ww_L(ni,l) * al * ww_L(nj,l)
                 
               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO
             
!              ENDIF

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_00_L_M


!------------------------------------------------------------------------------


SUBROUTINE qs_11_L_M (m0, jj, alpha, AA)
!=======================================

!  alpha << (Dw), (D_) >>   ===>   AA    incremental accumulation

   USE Gauss_points
   USE Gauss_points_L

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),                 INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: AA

   INTEGER :: mm, m, l, ni, nj, i, j, p
   REAL(KIND=8) :: al, x

   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G_L

         al = alpha * rj_L(l,m)

         DO ni = 1, n_w_L;  i = jj(ni, m)

            DO nj = 1, n_w_L;  j = jj(nj, m)

!              IF (j >= i) THEN
             
               x = al * SUM(dw_L(:,ni,l,m) * dw_L(:,nj,l,m))
              
               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO
           
!              ENDIF

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_11_L_M


END MODULE qs_L_M
