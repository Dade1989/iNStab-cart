!=========================
!=========================

PROGRAM  main
!
! Author: Jacopo Canton
! E-mail: jacopo.canton@mail.polimi.it
! Last revision: 30/8/2013
!
!=========================
!=========================

   USE ISO_C_BINDING
   USE loca_wrappers

   USE dynamic_structures
   USE prep_mesh_p1p2_sp
   USE sparse_matrix_profiles
   USE global_variables
   USE Dirichlet_Neumann
   USE Gauss_points
   USE Gauss_points_L
   USE start_sparse_kit
   USE qv_sp
   USE qc_sp_M
   USE par_solve_mumps
   USE vtk_plot
   USE restart_io
   USE cartesian_boundary_values
   USE read_input_files
   USE transient_growth
   USE dns_algorithms
   USE vorticity_stream
! Centre Manifold
   USE centreManifold
   USE cm_time_evolution
   USE qb_sp_m
   USE stressModule
! Harmonic Gain
   USE harmonicGain
! fem_boundary
   USE fem_boundary

!------------------------------------------------------------------------------


   IMPLICIT NONE

   REAL(KIND=8), PARAMETER :: zero = 0,  one = 1

! Test boundary connectivity ----------------------
   INTEGER , DIMENSION(:)   , ALLOCATABLE :: side_Id
   INTEGER , DIMENSION(:,:) , ALLOCATABLE :: bndEl
   INTEGER , DIMENSION(:)   , ALLOCATABLE :: iGlo
   INTEGER , DIMENSION(:,:) , ALLOCATABLE :: jjs5 , jj5Belt
   TYPE(CSR_MUMPS_Matrix) :: Mass_s , X_s , V_s
   INTEGER :: i1
! Stress computation ------------------------------
   REAL(KIND=8) , DIMENSION(:,:), ALLOCATABLE :: tauW
   REAL(KIND=8) , DIMENSION(:)  , ALLOCATABLE :: force
! -------------------------------------------------

   INTEGER      :: k, m
   REAL(KIND=8) :: dummy

!-------------END OF DECLARATIONS----------------------------------------------
!------------------------------------------------------------------------------
   WRITE(*,*) ''
   WRITE(*,*) ''
!------------------------------------------------------------------------------
!-------------INITIALIZE MPI---------------------------------------------------

   CALL MPI_INIT (mpiIerr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank, mpiIerr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nProc,  mpiIerr)
   IF ( mpiIerr .EQ. MPI_SUCCESS ) THEN
      IF ( myRank .EQ. 0 ) THEN
         WRITE(*,*) 'MPI correctly initialized'
         WRITE(*,*) 'number of processes:', nProc
      ENDIF
   ELSE
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** ERROR:                        ***'
      WRITE(*,*) '*** MPI uncorrectly initialized   ***'
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      STOP
   ENDIF

IF ( myRank == 0 ) THEN

!------------------------------------------------------------------------------
!-------------PROGRAM DATA INPUT-----------------------------------------------

   CALL read_program_data('program_data.in',  p_in)

!------------------------------------------------------------------------------
!-------------PREPARE P1/P2 GRID-----------------------------------------------

   CALL read_p1_gen_p2_sp (p_in%mesh_directory, p_in%mesh_name)
   
   number_of_sides = MAXVAL(sides)

!------------------------------------------------------------------------------
!-------------PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY---------------

   CALL read_and_apply_boundary_conditions('problem_data.in', k_d, rr, mms, jjs, &
                         sides,  Re, alphaV, js_D, zero_bvs_D, old_bvs_D, bvs_D, &
                                                          ms_2, ms_3, c_2, q_3)
   
!------------------------------------------------------------------------------
!------------ARRAY ALLOCATION--------------------------------------------------

   Nx = velCmpnnts*np + np_L
   ALLOCATE (uu(velCmpnnts, np),  u0(velCmpnnts, np))
   ALLOCATE (pp(np_L),            p0(np_L))
   ALLOCATE (xx(Nx),              x0(Nx))

!------------------------------------------------------------------------------
!------------MATRIX STRUCTURING ACCORDING TO CSR FORMAT------------------------

   CALL start_coupled_system (np, np_L, jj, js,  Jacobian)

   DESINGULARIZE = SIZE(ms_3) == 0  !  means Gamma_3 is void 

!------------------------------------------------------------------------------
!-------------SYMBOLIC FACTORIZATION OF MATRIX---------------------------------

   CALL par_mumps_master (INITIALIZATION, 1, Jacobian, 0)
   CALL par_mumps_master (SYMBO_FACTOR,   1, Jacobian, 0)
!
!------------------------------------------------------------------------------
!------------INITIAL GUESS EITHER FROM STOKES PROBLEM OR RESTART FILE----------
 
   IF ( p_in%read_restart_flag ) THEN

      CALL read_restart(x0, dummy, p_in%input_restart_file, LEN(trim(p_in%input_restart_file)))

      CALL extract(x0, u0, p0)

      ! Works only if Re is the only thing that changes from the restart
      ! file read and the new computation I want to make.
      ! No manual restart on parameters that change the boundary conditions
      ! has been implemented.
      ! Use LOCA instead.

   ELSE

      u0 = 0

      CALL compute_Stokes_initial_guess(np, mm, jj, jj_L, jjs, iis, js_D, bvs_D, &
                ms_2, ms_3, c_2, q_3, DESINGULARIZE, Re,   Jacobian, u0, p0, x0)

   END IF
   
   xx = x0


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------ANALYSIS----------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

   SELECT CASE ( p_in%method )
   CASE (1)
   !---------------------------
   ! Solve steady state
   !---------------------------

      ! my_null_ptr could be substituted with C_NULL_PTR?
      ite_num = nonlinear_solver_conwrap (xx, my_null_ptr, 1, Re, 0d0)

      CALL extract (xx,  uu, pp)

      IF ( p_in%write_restart_flag ) THEN
         ! WRITE RESTART FILE
         CALL write_restart(xx, Re, ite_num, p_in%nwtn_maxite, p_in%output_restart_file, LEN(trim(p_in%output_restart_file)))
      END IF

      ! WRITE QP RESTART FILE
      CALL write_QP_restart(xx, 'suiteSteadyState.QPrestart', 26)

      IF ( p_in%write_plots_flag ) THEN
         ! PLOT OUTPUT IN VTK FORMAT
         CALL vtk_plot_P2 (rr, jj, jj_L, uu, pp, trim(p_in%plot_directory) // 'steadyStateSolution.vtk')

         ! computation of vorticity and stream function
         ! set for the case of a single cylinder
         ALLOCATE (Dir_psi(number_of_sides))
         ALLOCATE (zz(np), psi(np))
         Dir_psi = (/.TRUE.,  &
                     .TRUE.,  &
                     .FALSE., &
                     .TRUE.,  &
                     .TRUE. /)
         CALL compute_vorticity_stream (mm, jj, js, uu, Dir_psi,  zz, psi)
         CALL vtk_plot_scalar_P2 (rr, jj,  zz, trim(p_in%plot_directory) // 'steadyStateVorticity.vtk')
         CALL vtk_plot_scalar_P2 (rr, jj, psi, trim(p_in%plot_directory) // 'steadyStateStream.vtk')
         DEALLOCATE(Dir_psi)
         DEALLOCATE(zz, psi)
      END IF

      IF ( p_in%write_BVS_flag ) THEN
         ! WRITE BOUNDARY VALUES TO FILE
         !CALL write_BVS (8, uu, rr, jjs, sides, 'steadyStateSolution')
         !CALL write_BVS (9, uu, rr, jjs, sides, 'steadyStateSolution')
      END IF


   CALL computeWallStress_real ( xx , rr , (/ 2 /) , tauW , force)
   WRITE(*,*)  
   WRITE(*,*) "Force = " , force
   WRITE(*,*)  

   CASE (2)
   !-------------------------------
   ! Continuation analisys
   !-------------------------------

      ! passdown structure
      pd%ldz      = Nx
      pd%x        = C_LOC(xx)
      pd%reynolds = Re
      pd%vRatio   = alphaV
      pd%mu       = 1d0
      pd%alpha    = 1d0
      pd%maxiter  = p_in%nwtn_maxite
      pd%tol      = p_in%nwtn_tol

      CALL do_loca(C_LOC(pd))

   CASE (3)
   !-------------------------------
   ! Eigenvalue computation on an
   ! already computed base flow
   !-------------------------------

      IF ( .NOT.p_in%read_restart_flag ) THEN
         WRITE(*,*)
         WRITE(*,*) '******************************************************'
         WRITE(*,*) '*** I suggest you don''t compute the eigenvalues    ***'
         WRITE(*,*) '*** of the Stokes initial guess. Load a previously ***'
         WRITE(*,*) '*** computed base flow.                            ***'
         WRITE(*,*) '******************************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      CALL compute_eigen(xx, 'SteadyState', 11, 0d0)

   CASE (4)
   !-------------------------------
   ! Structural sensitivity
   ! analysis on an already computed
   ! base flow, and already computed
   ! both direct and adjoint
   ! eigenvectors
   !-------------------------------

      CALL compute_structSens(Jacobian) ! Jacobian matrix only needed for its
                                        ! sparsity pattern
      
   CASE (5)
   !-------------------------------
   ! Transient growth computation on
   ! an already computed base flow
   !-------------------------------

      IF ( .NOT.p_in%read_restart_flag ) THEN
         WRITE(*,*)
         WRITE(*,*) '********************************************************'
         WRITE(*,*) '*** I suggest you don''t compute the transient growth ***'
         WRITE(*,*) '*** of the Stokes initial guess. Load a previously   ***'
         WRITE(*,*) '*** computed base flow.                              ***'
         WRITE(*,*) '********************************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      CALL compute_transientGrowth(x0, Jacobian, 'SteadyState')

   CASE (6)
   !-------------------------------
   ! DNS
   !-------------------------------
      CALL dns(x0)

   CASE (8)
   !-------------------------------
   ! Centre Manifold Reduction
   !-------------------------------
   ! Read the Base Flow, the parameter and the critical eigenspace
      CALL compute_centre_manifold()
   CASE (9)
      CALL evolution_on_centre_manifold()

   CASE (10)
   !-------------------------------
   ! Harmonic Gain
   !-------------------------------
      CALL compute_harmonic_gain(xx,Jacobian)

   CASE (11)
   !-------------------------------
   ! Stochastic Gain ONLY
   !-------------------------------
   ! read the Base Flow and the Optimal Forcing, Response and Gain
   !  and then compute the SENSITIVITY to:
   !  - Baseflow modifications
   !  - Steady Volume Forcing
   !  - Steady Wall Blowing/Suction
   !  - Control Cylinder

      WRITE(*,*) " TO BE WRITTEN YET ! "
      WRITE(*,*) " TO NOW SENSITIVITies ARE COMPUTED TOGETHER WITH "
      WRITE(*,*) " THE OPTIMAL RESPONSE IN CASE(10)                "
      STOP

   CASE DEFAULT

      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** Wrong parameter:              ***'
      WRITE(*,*) '*** p_in % method                 ***'
      WRITE(*,*) '*** set to: ', p_in%method
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      STOP

   END SELECT

!------------------------------------------------------------------------------
!-----------MPI FINALIZATION AND INFINITE LOOP FOR SLAVES----------------------

   CALL par_mumps_master (FINALIZATION, 1, Jacobian, 0)

   ELSE
      DO WHILE ( 1 > 0 )
         CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
         CALL par_mumps (parMumpsJob, matrID)
      ENDDO
   ENDIF

!==============================================================================
!==============================================================================
! END OF MAIN PROGRAM
!==============================================================================
!==============================================================================

CONTAINS

!------------------------------------------------------------------------------
! old subroutines

FUNCTION user_time () RESULT (u_t)

   ! NOT USED

   REAL :: u_t !! SINGLE PRECISION

   CALL CPU_TIME(u_t)

END FUNCTION user_time

!------------------------------------------------------------------------------

SUBROUTINE check_boundary_conditions (Dx, Dy, Norm, Tang)

   ! NOT USED

   IMPLICIT NONE

   LOGICAL, INTENT(IN) :: Dx, Dy, Norm, Tang 

   IF (Dx .AND. Dy) THEN ! fully vectorial Dirichlet condition

      IF (Norm .OR. Tang) THEN 
         
      write(*,*) '          NEWTON'
      write(*,*) '    SUB : check_boundary_conditions'
         WRITE (*,*) '--> With a fully vectorial Dirichlet condition'
         WRITE (*,*) '--> normal and/or tangential boundary conditions' 
         WRITE (*,*) '--> cannot be prescribed.  STOP.'
         WRITE (*,*) '--> Modify the input setting'
         WRITE (*,*) '--> Norm = .false.  and  Tang =.false..'
      write(*,*)
         STOP 

      ENDIF 

   ELSEIF (.NOT.Dx  .AND.  .NOT.Dy) THEN ! No vector Dirichlet condition

      IF (.NOT.Norm  .OR.  .NOT.Tang) THEN   
        write(*,*) '          NEWTON'
      write(*,*) '    SUB : check_boundary_conditions'
         WRITE (*,*) '--> SENSITE CONDITIONS REQUIRED'
         WRITE (*,*) '--> Norm and Tang must be both true.' 
         WRITE (*,*) '--> STOP.' 
         WRITE (*,*) '--> Modify the input setting'
         WRITE (*,*) '--> Norm = .true.  and  Tang =.true..'
      write(*,*)
         STOP 

      ENDIF

   ELSEIF ((Dx .AND. .NOT.Dy) .OR. (.NOT.Dx .AND. Dy)) THEN

      IF ((Norm .AND. Tang) .OR. (.NOT.Norm .AND. .NOT.Tang)) THEN

          write(*,*) '          NEWTON'
      write(*,*) '    SUB : check_boundary_conditions'          
         WRITE (*,*) '--> Only one beween Norm and Tang is permitted.'
         WRITE (*,*) '--> STOP.'
         WRITE (*,*) '--> The input must be modified.'
      write(*,*)
         STOP
         
      ENDIF

   ELSE 

       write(*,*) '          NEWTON'
     write(*,*) '    SUB : check_boundary_conditions'         
      WRITE (*,*) '--> NEVER HERE.  STOP.'  
     write(*,*) '--> hai creato un paradosso in barba'
     write(*,*) '--> alla logica.'

      STOP

   ENDIF 

END SUBROUTINE check_boundary_conditions

!-------------------------------------------------------------------------





!------------------------------------------------------------------------------
! new subroutines

SUBROUTINE read_and_apply_boundary_conditions(input_file, k_d, rr, mms, jjs, &
                     sides,  Re, alphaV, js_D, zero_bvs_D, old_bvs_D, bvs_D, &
                                                       ms_2, ms_3, c_2, q_3)

   IMPLICIT NONE
   ! input variables
   CHARACTER(*),                 INTENT(IN) :: input_file
   INTEGER,                      INTENT(IN) :: k_d
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: jjs
   INTEGER,      DIMENSION(:),   INTENT(IN) :: mms
   INTEGER,      DIMENSION(:),   INTENT(IN) :: sides
   ! output variables
   REAL(KIND=8)                        :: Re, alphaV
   TYPE(dyn_int_line),  DIMENSION(:)   :: js_D
   TYPE(dyn_real_line), DIMENSION(:)   :: zero_bvs_D, old_bvs_D, bvs_D
   INTEGER,      DIMENSION(:), POINTER :: ms_2, ms_3
   REAL(KIND=8), DIMENSION(:), POINTER :: c_2,  q_3
   ! local variables
   LOGICAL, ALLOCATABLE, DIMENSION(:)  :: Norm_u, Tang_u

! ----
   CHARACTER(LEN=1) , DIMENSION(:) , ALLOCATABLE :: id_Sides
! ----

   ! executable statements

   ALLOCATE ( Dir(velCmpnnts, number_of_sides),  & 
              Norm_u(number_of_sides),  Tang_u(number_of_sides) )
   ALLOCATE ( in_bvs_D(velCmpnnts, number_of_sides, 5) )



   ! (1)
   ! read input file
   OPEN (UNIT = 21, FILE = trim(input_file), FORM = 'formatted', STATUS = 'unknown')
   READ  (21,*)  Re, alphaV
   READ  (21,*) ! comment line
     
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> reading input-file ', trim(input_file), ' ...'
   WRITE(*,*) '    Re    = ', Re
   WRITE(*,*) '    alpha = ', alphaV
   WRITE(*,*) '    number of boundary segments: ', number_of_sides

! ----
  ALLOCATE(id_Sides(number_of_sides))
! ----

   DO m = 1, number_of_sides
  
      READ  (21,*)  Dir(1,m),  Dir(2,m),  Norm_u(m),  Tang_u(m) , id_Sides(m)
!      CALL check_boundary_conditions (Dir(1,m), Dir(2,m), Norm_u(m), Tang_u(m))
      WRITE (*,*) '    side m = ', m, ' ok'
  
   ENDDO

! ----
   write(*,*) "id_Sides = " , id_Sides
   CALL build_Id_sides_str ( id_Sides , sides_label_structure )
   write(*,*) "sides_label_structure : "
   IF (ALLOCATED(sides_label_structure%inflow )) THEN
     write(*,*) "     % inflow  = " , sides_label_structure%inflow
   ENDIF
   IF (ALLOCATED(sides_label_structure%wall   )) THEN
     write(*,*) "     % wall    = " , sides_label_structure%wall
   ENDIF
   IF (ALLOCATED(sides_label_structure%outflow)) THEN
     write(*,*) "     % outflow = " , sides_label_structure%outflow
   ENDIF
   IF (ALLOCATED(sides_label_structure%axis   )) THEN
     write(*,*) "     % axis    = " , sides_label_structure%axis
   ENDIF
!  build 'logical' matrix for Wal and DirInflow
   CALL build_Wal_DirInflow (velCmpnnts, number_of_sides , &
                             sides_label_structure , logicWa , logicIn )
! ----

   READ(21,*) ! jump one line
   ! boundary values for each side
   DO m = 1, number_of_sides

      READ (21,*) in_bvs_D(1,m,1), in_bvs_D(1,m,2), in_bvs_D(1,m,3), in_bvs_D(1,m,4), in_bvs_D(1,m,5), &
                  in_bvs_D(2,m,1), in_bvs_D(2,m,2), in_bvs_D(2,m,3), in_bvs_D(2,m,4), in_bvs_D(2,m,5)

   ENDDO
   
  
   CLOSE (21)
   WRITE(*,*) '--> finished reading file ', trim(input_file), '.'

   ! (2)
   ! generate Dirichlet nodes
   DO k = 1, velCmpnnts
   
      CALL Dirichlet_nodes_gen (jjs, sides, Dir(k,:),  js_D(k)%DIL)
      
      ALLOCATE (     bvs_D(k)%DRL(SIZE(js_D(k)%DIL)))
      ALLOCATE (zero_bvs_D(k)%DRL(SIZE(js_D(k)%DIL)))
      ALLOCATE ( old_bvs_D(k)%DRL(SIZE(js_D(k)%DIL)))

      zero_bvs_D(k)%DRL  = zero
  
   ENDDO

   ! (3)
   ! generate Dirichlet boundary values
   CALL gen_dirichlet_boundary_values (rr, sides, Dir, jjs, js_D, in_bvs_D, bvs_D)
   DO k = 1, velCmpnnts
      old_bvs_D(k)%DRL = bvs_D(k)%DRL
   ENDDO

   ! (4)
   ! generate Neumann boundary values
   ALLOCATE (c_2(SIZE(js)),  q_3(SIZE(js))) ! possibly SIZE = 0
   
   CALL Neumann_elements_gen (mms, sides, Norm_u,  ms_2)
   CALL Neumann_elements_gen (mms, sides, Tang_u,  ms_3)
 
   WRITE(*,*)
   IF (ANY(Norm_u)) THEN 
      
      WRITE (*,*) '    Gamma_2 is NOT void:', '  SIZE(ms_2) =', SIZE(ms_2) 
     
      c_2 = 0
   
   ENDIF
   
   IF (ANY(Tang_u)) THEN 
   
      WRITE (*,*) '    Gamma_3 is NOT void:', '  SIZE(ms_3) =', SIZE(ms_3)  
      
      q_3 = 0
   
   ENDIF

END SUBROUTINE read_and_apply_boundary_conditions

!------------------------------------------------------------------------------

SUBROUTINE compute_Stokes_initial_guess(np, mm, jj, jj_L, jjs, iis, js_D, bvs_D, &
                  ms_2, ms_3, c_2, q_3, DESINGULARIZE, Re,  Jacobian, u0, p0, x0)

   IMPLICIT NONE
   ! input variables
   INTEGER,                              INTENT(IN) :: np
   INTEGER, DIMENSION(:),                INTENT(IN) :: mm
   INTEGER, DIMENSION(:,:),              INTENT(IN) :: jj, jj_L, jjs, iis
   TYPE(dyn_int_line),  DIMENSION(:),    INTENT(IN) :: js_D
   TYPE(dyn_real_line), DIMENSION(:),    INTENT(IN) :: bvs_D
   INTEGER,      DIMENSION(:), POINTER,  INTENT(IN) :: ms_2, ms_3
   REAL(KIND=8), DIMENSION(:), POINTER,  INTENT(IN) :: c_2,  q_3
   LOGICAL,                              INTENT(IN) :: DESINGULARIZE
   REAL(KIND=8),                         INTENT(IN) :: Re
   ! output variables
   TYPE(CSR_MUMPS_Matrix)       :: Jacobian
   REAL(KIND=8), DIMENSION(:,:) :: u0
   REAL(KIND=8), DIMENSION(:)   :: p0
   REAL(KIND=8), DIMENSION(:)   :: x0
   ! local variables
   INTEGER :: Nx
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vv


   ! executable statements

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Computing Stokes initial guess'


   ALLOCATE ( vv(SIZE(u0,1), SIZE(u0,2)) )

   ! START OF STOKES INITIAL GUESS
   !==============================

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             Jacobian  <---  K_  +  V_ (weak)  
   !             Jacobian  <---  - V._   
   !-------------ONLY THE CONSTANT CONTRIBUTION---------------------- 

   WRITE(*,*) '*check*'
   WRITE(*,*) '    Re = ', Re
   CALL ComputeJacobianMatrix (np, mm, jj, jj_L, js_D, DESINGULARIZE, Jacobian, Re)

   CALL par_mumps_master (NUMER_FACTOR, 1, Jacobian, 0)
 
   !------------------------------------------------------------------
   !-------------COMPUTE THE RIGHT-HAND SIDE--------------------------
    
   vv = 0
   
   CALL qc_t0_sp_s (ms_2, jjs, iis,  c_2,  vv)  !  cumulative
   CALL qc_n0_sp_s (ms_3, jjs, iis, -q_3,  vv)  !  cumulative
   
   CALL collect (vv, 0*p0,  x0) ! here x0 is the RHS

   !------------------------------------------------------------------
   !-------------ENFORCING NONHOMOGENEOUS DIRICHLET BOUNDARY VALUES---

   CALL Dirichlet_c (np, js_D, bvs_D,  x0)
    
   Nx = SIZE(x0)
   IF (DESINGULARIZE) x0(Nx : Nx) = 0
                       
   !------------------------------------------------------------------
   !-------------DIRECT SOLUTION OF THE STOKES PROBLEM----------------

   CALL par_mumps_master (DIRECT_SOLUTION, 1, Jacobian, 0, x0)
   
   DEALLOCATE ( vv )


   CALL extract (x0,  u0, p0)

   WRITE (*,*) '    End of the Stokes initial guess'       
   WRITE (*,*)
   
!   CALL vtk_plot_P2 (rr, jj, jj_L, u0, p0, trim(p_in%plot_directory) // 'steadyStateSTOKES.vtk')

   ! END STOKES INITIAL GUESS
   !=========================

END SUBROUTINE compute_Stokes_initial_guess

!==============================================================================

END PROGRAM  main
