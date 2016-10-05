MODULE Gauss_points_L

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: k_d_L = 2,  n_w_L  = 3,  l_G_L  = 3,   &
                                             n_ws_L = 2,  l_Gs_L = 2
 
   REAL(KIND=8), DIMENSION(n_w_L,  l_G_L),           PUBLIC :: ww_L
   REAL(KIND=8), DIMENSION(n_ws_L, l_Gs_L),          PUBLIC :: wws_L

   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE, PUBLIC :: dw_L
   REAL(KIND=8), DIMENSION(:,    :, :), ALLOCATABLE, PUBLIC :: rnorms_L
   REAL(KIND=8), DIMENSION(:, :),       ALLOCATABLE, PUBLIC :: rj_L
   REAL(KIND=8), DIMENSION(:, :),       ALLOCATABLE, PUBLIC :: rjs_L
   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE, PUBLIC :: dwps_L !special!
   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE, PUBLIC :: dws_L  !SPECIAL!

  !REAL(KIND=8), DIMENSION(k_d_L,  n_w_L,  l_G_L,   me),  PUBLIC :: dw_L
  !REAL(KIND=8), DIMENSION(k_d_L,          l_Gs_L,  mes), PUBLIC :: rnorms_L
  !REAL(KIND=8), DIMENSION(l_G_L,   me),                  PUBLIC :: rj_L
  !REAL(KIND=8), DIMENSION(l_Gs_L,  mes),                 PUBLIC :: rjs_L

   PUBLIC Gauss_gen_L


CONTAINS


SUBROUTINE Gauss_gen_L (np, me, nps, mes, jj, jjs, rr)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: np, me, nps, mes

   INTEGER,      DIMENSION(n_w_L,  me),  INTENT(IN) :: jj
   INTEGER,      DIMENSION(n_ws_L, mes), INTENT(IN) :: jjs
   REAL(KIND=8), DIMENSION(k_d_L,  np),  INTENT(IN) :: rr

   REAL(KIND=8), DIMENSION(k_d_L, n_w_L,  l_G_L) :: dd
   REAL(KIND=8), DIMENSION( 1 ,  n_ws_L, l_Gs_L) :: dds
   REAL(KIND=8), DIMENSION(l_G_L)                :: pp
   REAL(KIND=8), DIMENSION(l_Gs_L)               :: pps

   REAL(KIND=8), DIMENSION(n_w_L)        :: r
   REAL(KIND=8), DIMENSION(n_ws_L)       :: rs
   REAL(KIND=8), DIMENSION(k_d_L, k_d_L) :: dr
   REAL(KIND=8), DIMENSION( 1 , k_d_L)   :: drs

   REAL(KIND=8) :: rj_Lac, rj_Lacs
   INTEGER      :: m, l, k, h, n,  ms, ls, dummy

   dummy = nps ! otherwise nps is not used 

   IF (ALLOCATED(dw_L)) THEN
    
      DEALLOCATE(dw_L, rnorms_L, rj_L, rjs_L, dwps_L)
  
   END IF

   ALLOCATE(dw_L(k_d_L, n_w_L, l_G_L,  me ))
   ALLOCATE(rnorms_L(k_d_L,  l_Gs_L, mes))
   ALLOCATE( rj_L(l_G_L,  me ))
   ALLOCATE(rjs_L(l_Gs_L, mes))
   ALLOCATE(dwps_L(1,  n_ws_L,  l_Gs_L,  mes))
   ALLOCATE(dws_L (1,  n_ws_L,  l_Gs_L,  mes))

!  evaluate and store the values of derivatives and of the
!  jacobian determinant at Gauss points of all volume elements

!  volume elements

   CALL element_2d (ww_L, dd, pp)

   DO m = 1, me

      DO l = 1, l_G_L

         DO k = 1, k_d_L
            r = rr(k, jj(:,m))
            DO h = 1, k_d_L
               dr(k, h) = SUM(r * dd(h,:,l))
            ENDDO
         ENDDO

         rj_Lac = dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1)

         DO n = 1, n_w_L
            dw_L(1, n, l, m)   &
               = (+ dd(1,n,l)*dr(2,2) - dd(2,n,l)*dr(2,1))/rj_Lac
            dw_L(2, n, l, m)   &
               = (- dd(1,n,l)*dr(1,2) + dd(2,n,l)*dr(1,1))/rj_Lac
         ENDDO

         rj_L(l, m) = rj_Lac * pp(l)

      ENDDO

   ENDDO


!  surface elements

   CALL element_1d (wws_L, dds, pps)

   DO ms = 1, mes

      DO ls = 1, l_Gs_L

         DO k = 1, k_d_L
            rs = rr(k, jjs(:,ms))
            drs(1, k) = SUM(rs * dds(1,:,ls))
         ENDDO

         rj_Lacs = SQRT( drs(1,1)**2 + drs(1,2)**2 )

         rnorms_L(1, ls, ms) = + drs(1,2)/rj_Lacs   ! outward normal
         rnorms_L(2, ls, ms) = - drs(1,1)/rj_Lacs   ! outward normal

         rjs_L(ls, ms) = rj_Lacs * pps(ls)

      ENDDO

   ENDDO


   DO ms = 1, mes  ! necessary only for evaluating gradient
                   ! tangential to the surface (ex COMMON Gauss_tan)

      DO ls = 1, l_Gs_L

          dwps_L(1, :, ls, ms) = dds(1, :, ls) * pps(ls)
           dws_L(1, :, ls, ms) = dds(1, :, ls)

      ENDDO

   ENDDO

!  write(*,*) '    MOD : Gauss_points_L'
!  write(*,*) '    SUB : Gauss_gen_L'
!  write(*,*) '    CAL : read_p1_gen_p2_sp <- prep_mesh_p1p2_sp'
   PRINT*,     '    end of gen_Gauss_L'
!    write(*,*)


   CONTAINS
   !=======

   SUBROUTINE element_2d (w, d, p)

!     triangular element with linear interpolation
!     and three Gauss integration points

!        w(n_w_L, l_G_L) : values of shape functions at Gauss points
!     d(2, n_w_L, l_G_L) : derivatives values of shape functions at Gauss points
!             p(l_G_L) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_w_L, l_G_L), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(2, n_w_L, l_G_L), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_G_L),           INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_G_L) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  six = 6

      REAL(KIND=8) :: f1, f2, f3, x, y
      f1(x, y) = one - x - y
      f2(x, y) = x
      f3(x, y) = y

      xx(1) = one/six;  xx(2) = two/three;  xx(3) = one/six
      yy(1) = one/six;  yy(2) = one/six;    yy(3) = two/three

      DO j = 1, l_G_L

            w(1, j) = f1(xx(j), yy(j))
         d(1, 1, j) = - one
         d(2, 1, j) = - one

            w(2, j) = f2(xx(j), yy(j))
         d(1, 2, j) = one
         d(2, 2, j) = zero

            w(3, j) = f3(xx(j), yy(j))
         d(1, 3, j) = zero
         d(2, 3, j) = one

               p(j) = one/six

      ENDDO

   END SUBROUTINE element_2d

!------------------------------------------------------------------------------

   SUBROUTINE element_1d (w, d, p)

!     one-dimensional element with linear interpolation
!     and two Gauss integration points

!        w(n_w_L, l_G_L) : values of shape functions at Gauss points
!     d(1, n_w_L, l_G_L) : derivatives values of shape functions at Gauss points
!             p(l_G_L) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_ws_L, l_Gs_L), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(1, n_ws_L, l_Gs_L), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_Gs_L),            INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_Gs_L) :: xx
      INTEGER :: j

      REAL(KIND=8) :: one = 1,  two = 2,  three = 3

      REAL(KIND=8) :: f1, f2, x
      f1(x) = (one - x)/two
      f2(x) = (x + one)/two

      xx(1) = - one/SQRT(three)
      xx(2) = + one/SQRT(three)

      DO j = 1, l_Gs_L

            w(1, j) = f1(xx(j))
         d(1, 1, j) = - one/two

            w(2, j) = f2(xx(j))
         d(1, 2, j) = + one/two

               p(j) = one

      ENDDO

   END SUBROUTINE element_1d


END SUBROUTINE Gauss_gen_L

END MODULE Gauss_points_L
