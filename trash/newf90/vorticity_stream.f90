MODULE vorticity_stream

   USE global_variables
   USE prep_mesh_p1p2_sp ! for some global variables as jj
   USE Dirichlet_Neumann ! for Dirichlet_nodes_gen subroutine
   USE start_sparse_kit  ! for start_matrix_2d_p2
   USE qs_sp
   USE qs_sp_M
   USE par_solve_mumps

   IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE  compute_vorticity_stream (mm, jj, js, uu, Dir_psi,  zz, psi)

!  Compute the vorticity field  zz  and the stream function  psi
!  corresponding to the 2D solenoidal velocity field  uu

!  If the 2D domain is multiply connected, i.e., when there is one
!  or more immersed bodies in the flow, the program determines the 
!  constante values of the stream function by evaluating the 
!  flux through a path connecting one point on the boundary with
!  the immersed bodies.

!  rr, sides, pp:  GLOBAL VARIABLES

!  WARNING!! This program must be adapted to match the number of 
!            internal boundaries.
!            The program lines to be modified are
!            written between the marker: !***

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),     INTENT(IN) :: mm, js
   INTEGER,      DIMENSION(:,:),   INTENT(IN) :: jj
   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN) :: uu
   LOGICAL,      DIMENSION(:),     INTENT(IN) :: Dir_psi
   REAL(KIND=8), DIMENSION(:),     INTENT(OUT):: zz, psi

   TYPE(CSR_MUMPS_Matrix) :: MK  !  Mass & Stiffness matrix

   INTEGER,      DIMENSION(:), POINTER     :: js_psi_D
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: as_psi_D

   TYPE(dyn_int_line), DIMENSION(:),   ALLOCATABLE :: js_a
   REAL(KIND=8),       DIMENSION(:),   ALLOCATABLE :: psi_a
   REAL(KIND=8),       DIMENSION(:,:), ALLOCATABLE :: rc_a ! center of airfoil
   LOGICAL,            DIMENSION(SIZE(Dir_psi))    :: profile_l

   REAL(KIND=8),       DIMENSION(:),   ALLOCATABLE :: dist_vect

   REAL(KIND=8) :: psi_0, psi_l
   INTEGER      :: loops, np_L,  l, n, ns, n_0, el_0, j, k
   REAL(KIND=8) :: x_min

   INTEGER, DIMENSION(1) :: dummy


!------------------------------------------------------------------------------
!-------------VORTICITY COMPUTATION--------------------------------------------

   WRITE (*,*) ' Structuring of the matrix for vorticity and stream problems'

   CALL start_matrix_2d_p2 (SIZE(uu,2), jj, js,  MK)

   CALL par_mumps_master (INITIALIZATION, 8, MK, 0)
   CALL par_mumps_master (SYMBO_FACTOR,   8, MK, 0)
   
   WRITE (*,*) ' symbolic factorization of MM or KK matrices '
  
   ALLOCATE (MK%e(SIZE(MK%j)))

   ! right hand side for both the vorticity equation and
   ! the stream function Poisson equation

   CALL qs_01_sp_c (mm, jj, uu,  zz)  !  zz = (w, k.Rot u)

   psi = zz  !  save the weak vorticity also for the psi equation

   ! determination of the vorticity nodal values

   MK%e = 0
   
   CALL qs_00_sp_M (mm, jj, 1.d0,  MK)

   CALL par_mumps_master (NUMER_FACTOR, 8, MK, 0)
   
   WRITE (*,*) ' numerical factorization of  MM zz = rot u  problem '
  
   CALL par_mumps_master (DIRECT_SOLUTION, 8, MK, 0, zz)
   
   WRITE (*,*) ' direct solution of  MM zz = rot u  problem '
  
   WRITE (*,*) 'Vorticity field computed'


!------------------------------------------------------------------------------
!-------------STREAM FUNCTION COMPUTATION--------------------------------------

   CALL Dirichlet_nodes_gen (jjs, sides, Dir_psi,  js_psi_D)

   ALLOCATE (as_psi_D(SIZE(js_psi_D)))


   MK%e = 0

   CALL qs_11_sp_M (mm, jj, 1.d0,  MK)

   CALL Dirichlet_M (js_psi_D,  MK)

   CALL par_mumps_master (NUMER_FACTOR, 8, MK, 0)
  
   WRITE (*,*) ' numerical factorization of  KK psi = zz  problem '
  
   CALL stream_bound_values (js_psi_D,  as_psi_D,  rr)


   CALL Dirichlet (js_psi_D, as_psi_D,  psi)


!------------------------------------------------------------------------------
!-------------PROCEDURE IN CASE OF INTERNAL BOUNDARIES-------------------------

   ! Use of Euler relation for counting the number of
   ! possible internal loops
   !
   ! loops = np - np_L - (me + np_L - 1)

   np_L = SIZE(pp)  !  pp  is a GLOBAL VARIABLE

   loops = SIZE(uu,2) - np_L - (SIZE(jj,2) + np_L - 1)

   WRITE (*,*) 'loops = ', loops 

   IF (loops > 0) THEN

      ALLOCATE (js_a(loops),  psi_a(loops),  rc_a(1,loops))

      
!***
      rc_a(:, 1) = (/0.0, 0.0/)      ! point inside CYLINDER
!***


      ! Find a node n_0 on the external boundary (on the left vertical side)
      ! at which the value of psi is prescribed
    
!      DO ns = 1, SIZE(js_psi_D);  n = js_psi_D(ns)
!         IF (rr(1,n) < -4.999) THEN;  n_0 = n;  EXIT;  ENDIF
!      ENDDO

!      n_0 = MINLOC(rr(1, js_psi_D))

      x_min = 1.0d40
      DO ns = 1, SIZE(js_psi_D); j = js_psi_D(ns)
        IF ( rr(1, j) < x_min ) THEN
          x_min = rr(1, j)
          n_0 = j
        ENDIF
      ENDDO

      psi_0 = rr(2, n_0)
      
!      WRITE(*,*)
!      WRITE(*,*) '--> Nodo in cui viene impostato psi_0: '
!      WRITE(*,*) '--> coord_x = ', rr(1,n_0)
!      WRITE(*,*) '--> coord_y = ', rr(2,n_0)
!      WRITE(*,*) '--> n_0     = ', n_0
!      WRITE(*,*) '--> psi_0   = ', psi_0

      ! Find an element el_0 containing the node n_0

      DO n = 1, SIZE(jj, 2)
         IF (PRODUCT(SPREAD(n_0, 1, SIZE(jj, 1)) - jj(:, n)) == 0) THEN
         el_0 = n;  EXIT
         ENDIF
      ENDDO


      !OPEN (UNIT = 27, FILE = 'path', FORM = 'formatted', STATUS = 'unknown')


      ! COMPUTATION OF THE UNKNOW VALUES OF PSI ON THE PROFILES

      DO l = 1, loops

         ! Assign the segment belonging to the l-th profile to TRUE
         ! in order to recognize all its nodes, for the purpose of
         ! terminating the path on the target profile

         ! This segment depends on how many sides constitute
         ! the various profiles

         profile_l = .FALSE.

         SELECT CASE (l)

!***    
            CASE (1)   ! the boundary of each internal plate
                       ! is made of four sides (Grassi's 
                       ! experiment with sharp edged plates)
            
                k = 5

                profile_l(k) = .TRUE.
!***

         END SELECT


         CALL Dirichlet_nodes_gen (jjs, sides, profile_l,  js_a(l)%DIL)

         IF (l /= 1) THEN

            ALLOCATE(dist_vect(SIZE(js_a(l-1)%DIL)))

            DO ns = 1, SIZE(js_a(l-1)%DIL)
               dist_vect(ns) = SUM((rr(:, js_a(l - 1)%DIL(ns)) - rc_a(:, l))**2)
            ENDDO

            dummy = MINLOC(dist_vect)


            DEALLOCATE(dist_vect)


            n_0 = js_a(l-1)%DIL(dummy(1))   ! n_0 is the closest node
                                            ! (on the (l-1)th boundary) 
                                            ! to the center of the l-th plate

            DO n = 1, SIZE(jj, 2)

               IF (PRODUCT(jj(:, n) - SPREAD(n_0, 1, SIZE(jj, 1))) == 0) THEN

                  ALLOCATE(dist_vect(3))

                  DO ns =1, 3
                     dist_vect(ns) = SUM((rr(:, jj(ns, n)) - rc_a(:, l))**2)
                  ENDDO

                  dummy = MINLOC(dist_vect)

                  DEALLOCATE(dist_vect)

                  IF (n_0 /= jj(dummy(1), n)) THEN;  el_0 = n;  EXIT;  ENDIF             
                  
                  ! el_0 is the closest element to the center of
                  ! the l-th plate among the ones containing n_0

               ENDIF

            ENDDO

         ENDIF


         CALL psi_on_airfoils (psi_0, n_0, el_0, rr, uu, jj, js_a(l)%DIL, neigh, &
                               rc_a(:, l), psi_l)

         psi_a(l) = psi_l;   psi_0 = psi_l

WRITE(*,*) 'l = ', l, 'psi_a(l) = ', psi_a(l)

         CALL Dirichlet (js_a(l)%DIL, SPREAD(psi_a(l), 1, SIZE(js_a(l)%DIL)),  psi)

      ENDDO

   ENDIF

!-------------END OF PROCEDURE FOR INTERNAL BOUNDARIES-------------------------
!------------------------------------------------------------------------------

   CALL par_mumps_master (DIRECT_SOLUTION, 8, MK, 0, psi)

   IF ( loops > 0 ) THEN
      psi = psi - psi_a(1) ! only needed in case of internal boundaries!!!
   END IF
   
   WRITE (*,*) ' direct solution of  KK psi = zz  problem '
 
   WRITE (*,*) 'Stream function computed'


   DEALLOCATE (MK%i, MK%i_mumps, MK%j, MK%e)
   CALL par_mumps_master (DEALLOCATION, 8, MK, 0)


END SUBROUTINE  compute_vorticity_stream

!-----------------------------------------------------------------------------

SUBROUTINE  psi_on_airfoils (p_i, j_i, e_i, rr, uu, jj, js_a, neigh, r_a,  p_f)

   IMPLICIT NONE

   REAL(KIND=8),                 INTENT(IN)  :: p_i
   INTEGER,                      INTENT(IN)  :: j_i, e_i
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: rr, uu
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj, neigh
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: js_a
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: r_a
   REAL(KIND=8),                 INTENT(OUT) :: p_f

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: rs_vect

   INTEGER,      DIMENSION(1000) :: path
   REAL(KIND=8), DIMENSION(2)    :: dr, norm
   INTEGER,      DIMENSION(50)   :: bb_j, bb_e

   REAL(KIND=8) :: flux
   INTEGER      :: l_p, l, j, j1, j2, j0, k, m, n, n1, n2, e, e_loc, counter

   INTEGER, DIMENSION(1) :: dummy ! JUST FOR STUPID PROGRAMMERS


   ! determination of the the first step from the boundary 
   ! the final point of the first step must not be on the boundary 
 
   l = 1
   j1 = j_i

   dummy = MINLOC(neigh(:, e_i))
   k = dummy(1)

   IF (neigh(k, e_i) == 0) THEN
      j2 = jj(k, e_i)
   ELSE
      j2 = jj(1, e_i)
      IF (j2 == j1)  j2 = jj(2, e_i)
   ENDIF



   DO n = 1, 3
      IF (jj(n, e_i) == j1) THEN;  n1 = n;  EXIT;  ENDIF
   ENDDO

   DO n = 1, 3
      IF (jj(n, e_i) == j2) THEN;  n2 = n;  EXIT;  ENDIF
   ENDDO



   IF (n1 + n2 == 3) j0 = jj(6, e_i)
   IF (n1 + n2 == 4) j0 = jj(5, e_i)
   IF (n1 + n2 == 5) j0 = jj(4, e_i)

   path(l)     = j1
   path(l + 1) = j0
   path(l + 2) = j2

   l = l + 2

   !-----------------------------------------------------------------------

   e = e_i

   endless: DO WHILE (.TRUE.)

      j = path(l)

      IF (l > 1e6) THEN
         WRITE(*, *) 'the path is longer than 1000 edges'
         STOP
      ENDIF

      DO n = 1, SIZE(js_a)
         IF (j == js_a(n)) EXIT endless
      ENDDO


      ! creation of the list of elements containing the node j

      e_loc = e
      counter = 1

      DO WHILE (.TRUE.)

         ! position of the bubble center in the local element

         DO n1 = 1, 3
            IF (jj(n1, e_loc) == j) THEN;  k = n1;  EXIT;  ENDIF
         ENDDO                                                                    
                                             ! k = 1  -->  m = 3
                                             ! k = 2  -->  m = 1
         m = MODULO(k, 3) + 2                ! k = 3  -->  m = 2

         IF (m == 4) m = 1

         bb_j(counter) = jj(m, e_loc)
         bb_e(counter) = e_loc

         e_loc = neigh(MODULO(k, 3) + 1, e_loc)

         IF (e_loc == e) EXIT

         counter = counter + 1

      ENDDO


      ! vector of the distances from the bubble points to rc_a

      ALLOCATE (rs_vect(counter))

      ! rs_vect = SUM((rr(:, bb_j(1:counter)) - SPREAD(r_a, 2, counter))**2)

      DO n = 1, counter
         rs_vect(n) = SUM((rr(:, bb_j(n)) - r_a)**2)
      ENDDO

      ! element update
      
      dummy = MINLOC(rs_vect)

      DEALLOCATE (rs_vect)

      e = bb_e(dummy(1))

      DO n = 1, 3
         IF (jj(n, e) == j) THEN;  n1 = n;  EXIT;  ENDIF
      ENDDO

      DO n = 1, 3
         IF (jj(n, e) == bb_j(dummy(1))) THEN;  n2 = n;  EXIT;  ENDIF
      ENDDO

      IF (n1 + n2 == 3) j0 = jj(6, e)
      IF (n1 + n2 == 4) j0 = jj(5, e)
      IF (n1 + n2 == 5) j0 = jj(4, e)

      l = l + 1;  path(l) = j0
      l = l + 1;  path(l) = bb_j(dummy(1))

   ENDDO endless

   l_p = l  ! lenght of the path


   flux = 0

   DO l = 1,  l_p - 2,  2

      j1 = path(l);  j0 = path(l+1);  j2 = path(l+2)

      dr = rr(:,j2) - rr(:,j1)
      norm(1) = dr(2);   norm(2) = - dr(1)

      flux = flux + SUM( (uu(:,j1) + 4*uu(:,j0) + uu(:,j2)) * norm )

   ENDDO

   p_f = p_i + flux / 6


END SUBROUTINE  psi_on_airfoils

!------------------------------------------------------------------------------

SUBROUTINE  stream_bound_values (js_s,  bvs_s,  rr)

!  This program defines boundary values for stream function 
!  for the CYLINDER problem

!  bvs_s  contains the (Dirichlet) boundary values
!  of the stream function psi. 

   IMPLICIT NONE
  
   INTEGER,      DIMENSION(:),   INTENT(IN)  ::  js_s
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: bvs_s
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: rr

   REAL(KIND=8) :: a, b, c, d, R

   
!  The boundary condition on the internal boundaries can be either 
!  a homogeneous Neumann condition, in which case no boundary 
!  values of psi is required, or
!  a NONhomgeneous Dirichlet condition (with a different constant
!  but UNKNOWN value on each internal boundaries), in which case
!  the stream value is determined by integrating the mass flowing 
!  between the internal boundaries. 
!   
!  The boundary condition on the outlet boundary (at  x = b)
!  is a (possibly NONhomogeneous Neumann) condition and therefore
!  no boundary value for psi is required there.
!
!  The condition on the inlet boundary (at  x = a)
!  and on the bottom and top sides (at y = c  and  y = d) 
!  of the external boundary are instead of Dirichlet type, 
!  with boundary values in general different from zero.  
 
!  The boundary values of the stream function are obtained
!  from the stream function of constant profile of the horizontal
!  velocity component.  This gives a linear dependance of psi 
!  on the vertical cartesian coordinate y, precisely:
 

a = MINVAL(rr(1,:))
b = MAXVAL(rr(1,:))
c = MINVAL(rr(2,:))
d = MAXVAL(rr(2,:))

R = 0.5d0 ! cylinder radius
   
   
   WHERE ( rr(1,js_s)**2 + rr(2,js_s)**2 > 1.001*R**2 )
      
      bvs_s  =  1d0 * rr(2,js_s)

   END WHERE
   
  
END SUBROUTINE  stream_bound_values

!==============================================================================

END MODULE vorticity_stream
