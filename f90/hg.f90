!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
   USE dynamic_structures
   USE sparse_matrix_profiles
   USE sparse_matrix_operations
   USE global_variables
   USE miscellaneous_subroutines
   USE prep_mesh_p1p2_sp ! for some global variables as jj
   USE Dirichlet_Neumann ! for Dirichlet_nodes_gen subroutine
   USE qc_sp_M
   USE qv_sp
   USE qs_l_m
   USE qs_l_sp
   USE par_solve_mumps
   USE vtk_plot
   USE eigensolve ! needed if one wants to use an eigenvector as initial guess
   USE dns_algorithms
   USE axisym_boundary_values
   
   USE extendedModule

   IMPLICIT NONE
   ! variables "global" to this module
   TYPE(CSR_MUMPS_Complex_Matrix)               :: Wd, Zd, Wa, Za, W0, Z0, MassV
   TYPE(CSR_MUMPS_Complex_Matrix)               :: Gomega , Qf , Q , Bf , Cu_cmplx , CuGain , BBf , Bfff
! ---------------------------------------------------------------- Matrices for HG for Input Forcing -----
   TYPE(CSR_MUMPS_Complex_Matrix)               :: Qff_In , HQff , HQff_T , Qf_in , HtimesQ
! ---------------------------------------------------------------- Matrices for control sensitivity ------
   TYPE(CSR_MUMPS_COMPLEX_MATRIX)               :: eLNS , eeLNS , mLagr1 , mLagr , Du , DuEss
   TYPE(CSR_MUMPS_COMPLEX_MATRIX)               :: qf1 , qf3
! --------------------------------------------------------------------------------------------------------
   COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: uu_tg, u0_tg
   COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:)   :: pp_tg, p0_tg
   COMPLEX(KIND=8), DIMENSION(:), POINTER       :: xx_tg, x0_tg
   COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: u0opt_cmplx
   LOGICAL, ALLOCATABLE, DIMENSION(:,:)         :: Dir_tg
   TYPE(dyn_int_line), DIMENSION(velCmpnnts)    :: js_D_tg
   LOGICAL                                      :: DESINGULARIZE_tg
   TYPE(dyn_real_line), DIMENSION(velCmpnnts)   :: zero_bvs_D_tg
   REAL(KIND=8) :: dtE, dtH

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE computeHarmonicGain(x_vec, Lns, filenm)

   IMPLICIT NONE

   ! input variables
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x_vec  ! computed base flow
   TYPE(CSR_MUMPS_Matrix)                 :: Lns    ! Linearized Navier--Stokes operator (Jacobian matrix)
   CHARACTER(*),               INTENT(IN) :: filenm
   ! "output" variables
   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: singularValue  ! these two have a dimension more than necessary
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: singularVector ! because the algorithm for the SVD is "general"
   COMPLEX(KIND=8), DIMENSION(np)               :: Ek0  , Ek      ! kinetic energy fields
   COMPLEX(KIND=8)                              :: Ek0_s, Ek_s    ! kinetic energy
   ! local variables
   INTEGER      :: nsv
   COMPLEX(KIND=8) :: sigma
   INTEGER      :: nSteps
   INTEGER      :: i, k, ii ,i2 , i3 , i1
   LOGICAL      :: existFlag
   COMPLEX(KIND=8) :: resLinfty
   INTEGER      :: iteNum
   COMPLEX(KIND=8), DIMENSION(velCmpnnts,np) :: uInitGuess

   CHARACTER(LEN=128) :: restart_name
   REAL(KIND=8) :: dummy

   REAL(KIND=8) :: dataIdentifier = 313d0 ! Donald Duck's plate number!
   INTEGER      :: tausNumber = 1
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: taus
   CHARACTER(LEN=128) :: tauName ! used to insert tau in file names
   INTEGER            :: tauNameTemp

   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmpVector
   
   REAL(KIND=8), DIMENSION(:),  ALLOCATABLE :: omega
   COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: y1 , y2 , y3 , y4 , y5 , y6 , y7
!   COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: rhsr1 , rhsr2
   COMPLEX(KIND=8) :: residuo1 , residuo2
 ! ZNAUPD ----------------------------------------------------------------------
   INTEGER :: ido , ncv , nev , lworkl , info , n , ldv
   CHARACTER(LEN=1) :: bmat
   CHARACTER(LEN=2) :: which
   INTEGER, DIMENSION(11) :: iparam
   INTEGER, DIMENSION(14) :: ipntr
   REAL(KIND=8) :: tol
   COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: workd , resid , workl 
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: v
   REAL(KIND=8)   , DIMENSION(:)  , ALLOCATABLE :: rwork
 ! ZNEUPD parameters -----------------------------------------------------------
   LOGICAL, DIMENSION(:), ALLOCATABLE :: lselect
   COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: d
!   COMPLEX(KIND=8) :: sigma
   COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: workev
   INTEGER :: ierr
   REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: rd
   LOGICAL          :: rvec   
   INTEGER :: maxit = 100
   
   COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: singularValues            !  G^2_k
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: singularVectors           !  [ ... |  f_k  | ... ]
   COMPLEX(KIND=8), DIMENSION(:  ), ALLOCATABLE :: singularVectorsP          ! (u_k p_k)
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: optimalResponse           !  [ ... |  u_k  | ... ]
   COMPLEX(KIND=8), DIMENSION(:  ), ALLOCATABLE :: adjointsingularVectorsP   ! (u^+_k; p^+_k)               [ ... | (u^+_k; p^+_k) | ... ]
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: adjointOptimalResponse    !  [ ... | u^+_k | ... ]

   COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: baseFlowSens_cmplx
   REAL   (KIND=8), DIMENSION(:)  , ALLOCATABLE :: baseFlowSensVct   
   REAL   (KIND=8), DIMENSION(:,:), ALLOCATABLE :: baseFlowSens
   
   INTEGER , DIMENSION(:) , ALLOCATABLE :: labelSides_Wall    , mmWall
   INTEGER , DIMENSION(:) , ALLOCATABLE :: labelSides_Neumann , mmNeumann
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: optimalResponseK
   COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: fu1 , fu2
   COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: RHSadjU , adjointBaseFlow
   INTEGER , DIMENSION(:)   , ALLOCATABLE :: MMSWall
   INTEGER , DIMENSION(:,:) , ALLOCATABLE :: JJSWall , JJSWallFromOne
   INTEGER , DIMENSION(:)   , ALLOCATABLE :: MMSInflow
   INTEGER , DIMENSION(:,:) , ALLOCATABLE :: JJSInflow , JJSInflowFromOne
   
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: adjointBaseU
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: wallBlowSuctSens
   
! Stochastic gain
   REAL   (KIND=8), DIMENSION(:)  , ALLOCATABLE :: sumGains
   REAL   (KIND=8), DIMENSION(:,:), ALLOCATABLE :: sumBaseFlowSens
   REAL   (KIND=8), DIMENSION(:,:), ALLOCATABLE :: sumVolForceSens
   REAL   (KIND=8), DIMENSION(:,:), ALLOCATABLE :: sumWallBlowSuctSens 
   REAL   (KIND=8)                              :: stocGain
   REAL   (KIND=8), DIMENSION(:)  , ALLOCATABLE :: stocBaseFlowSens
   REAL   (KIND=8), DIMENSION(:)  , ALLOCATABLE :: stocVolForceSens
   REAL   (KIND=8), DIMENSION(:)  , ALLOCATABLE :: stocWallBlowSuctSens 
  
   
   INTEGER , DIMENSION(:) , ALLOCATABLE :: js_InflowAll
   
!! ONCE IT WORKS   ! the 'adjoint'//... arrays are thought for multiple simultaneous solution of the adjoint linear(ized) system
!!                 ! forced by the optimal response array
   INTEGER :: nconv
   real(kind=8) :: realpart
   
   INTEGER , DIMENSION(:) ,  ALLOCATABLE :: iRowQf , iColQf
   
   CHARACTER(LEN=1024) :: filename
   CHARACTER(LEN=1024) :: format_string
   CHARACTER(LEN=1024) :: omegaString
   CHARACTER(LEN=1024) :: numGainString
   CHARACTER(LEN=1024) :: omegaFolder
   
   INTEGER :: numHGplot
   INTEGER :: sizeOmega
   INTEGER :: aERR
   INTEGER :: bothSensId
   
   COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: uuInlet
   CHARACTER(LEN=1024) :: NUMeIG


!----------------------------------------------------
   OPEN(UNIT=20, FILE='./HarmonicGain/omega.in')
   READ(20,*) sizeOmega
   ALLOCATE(omega(sizeOmega))
   DO i2 = 1 , sizeOmega
      READ(20,*) omega(i2)
   END DO      
   CLOSE(20)
   write(*,*) omega
!----------------------------------------------------
   WRITE(*,*) ''
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to computeHarmonicGain'
   WRITE(*,*) ''
!----------------------------------------------------

   ALLOCATE (uu_tg(velCmpnnts, np), u0_tg(velCmpnnts, np))
   ALLOCATE (pp_tg(np_L),           p0_tg(np_L))
   ALLOCATE (xx_tg(Nx),             x0_tg(Nx))

!--------------------------------------------------------
! PREPARE MATRICES FOR THE HARMONIC GAIN COMPUTATION

   ! (1)
   ! prepare boundary conditions
   !

   ALLOCATE ( Dir_tg(velCmpnnts, number_of_sides) )
   IF ( p_in%tranGrowth_BC == 1 ) THEN
      ! homogeneous Dirichlet on every border
      WRITE(*,*) '    boundary conditions: zero velocity on every border'
      Dir_tg = .TRUE.
      DESINGULARIZE_tg = .TRUE.
   ELSEIF ( p_in%tranGrowth_BC == 2 ) THEN
      ! same BCs as for the base flow but homogeneous
      ! Dirichlet where the base flow has nonhomogeneous Dirichlet
      WRITE(*,*) '    boundary conditions: same as base flow'
      Dir_tg = Dir
      DESINGULARIZE_tg = DESINGULARIZE
   ELSE
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** Wrong parameter:              ***'
      WRITE(*,*) '*** p_in % tranGrowth_BC          ***'
      WRITE(*,*) '*** set to: ', p_in%tranGrowth_BC
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF
   
   
   
   !---------------------------------------------------------------------------------------------------------
   ! Harmonic Gain for forcing as an inlet condition
   IF ( TRIM(p_in%HG_vol_in) == 'in' ) THEN
! Old way to compute Qf_In : extremely inefficient !!!      
!      ALLOCATE( Qff_In%i      (SIZE(Lns%i))       ); Qff_In%i       = Lns%i
!      ALLOCATE( Qff_In%i_mumps(SIZE(Lns%i_mumps)) ); Qff_In%i_mumps = Lns%i_mumps
!      ALLOCATE( Qff_In%j      (SIZE(Lns%j))       ); Qff_In%j       = Lns%j
!      ALLOCATE( Qff_In%e      (SIZE(Lns%e))       ); Qff_In%e       = 0.0d0
!      
!      WRITE(*,*) 'Building matrix q00_Inflow ... '
!      CALL qc_0y0_zero_sp_boundary_cmplx (mms_In, jj, jjs, 1.0d0,  Qff_In) ! mms_In contains the boundary elements on the Inflow sides
!      WRITE(*,*) '   Matrix q0y0_Inflow built !'
!      
!      
!      ALLOCATE(iRowQf(3*SIZE(nodes_Inflow_Raw)),iColQf(3*SIZE(nodes_Inflow_Raw)))
!      DO i2 = 1 , 3
!        iRowQf(1+(i2-1)*SIZE(nodes_Inflow_Raw):i2*SIZE(nodes_Inflow_Raw)) = nodes_Inflow_Raw + np*(i2-1)     ! -1 if CSR starts counting
!      END DO                                                                                      !  from 0
!      
!      iColQf = iRowQf
!!      WRITE(*,*) iColQf
!      
!      
!      WRITE(*,*) 'Extracting Qf_In from Qff_In ... '
!      CALL extractMatrix ( Qff_In , iRowQf , iColQf , Qf_In )
!      DEALLOCATE(Qff_In%e,Qff_In%j,Qff_In%i,Qff_in%i_mumps)
!      WRITE(*,*) '   Qf_In built !'
!      CALL writeSparseMatrixToFile ( Qf_In )
      
      
  ! build the matrix Qf_in (efficient way ? )!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL boundaryElementsOfSides ( mms , sides   , sides_label_structure%inflow , mmsInflow )
     CALL fromJJStoJJSw ( jjs , mmsInflow , JJSInflow , JJSInflowFromOne  )
     
     CALL start_matrix_P2SP2S_cmplx     ( JJSInflow , qf1 )
     CALL start_matrix_3d_p2SP2S_cmplx  ( JJSInflowFromOne , qf1 , Qf_in )
     write(*,*) "qf_in%i" , qf_in%i
     write(*,*) "qf_in%j" , qf_in%j
     write(*,*) "qf_in%i_mumps" , qf_in%i_mumps
     CALL qc_0y0_WallOnly_sp_boundary_cmplx  ( mmsInflow , jj , jjs , 1.0d0 , Qf_in )
  !   CALL writeSparseMatrixToFile ( Qf_In )
      
      CALL par_mumps_master (INITIALIZATION, 14, Qf_In, 0)
      CALL par_mumps_master (SYMBO_FACTOR,   14, Qf_In, 0)
      write(*,*) 'Inlet forcing: starting numerical factorization of Qf ...'
      CALL par_mumps_master (NUMER_FACTOR,   14, Qf_In, 0)
      write(*,*) '   numerical factorization of Qf: done !'
      
    ! build the matrix BBf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      CALL hgInletVelocityB ( np , 3 , nodes_Inflow , nodes_Inflow_Raw , BBf )
      
      write(*,*) 'size(BBf%e) , minval , maxval = ' , size(BBf%e), minval(dble(BBf%e)), maxval(dble(BBf%e))
      write(*,*) 'size(BBf%i) , minval , maxval = ' , size(BBf%i), minval(     BBf%i ), maxval(     BBf%i )
      write(*,*) 'size(BBf%j) , minval , maxval = ' , size(BBf%j), minval(     BBf%j ), maxval(     BBf%j )
      
   END IF
   
!   STOP
!---------------------------------------------------------------------------------------------------------
   
!!!! ????????????????????????????????????????????????????????????????????????????????????????
   DO k = 1, velCmpnnts
      CALL Dirichlet_nodes_gen (jjs, sides, Dir_tg(k,:) .AND. .NOT.Axis,  js_D_tg(k)%DIL)
      ALLOCATE (zero_bvs_D_tg(k)%DRL(SIZE(js_D_tg(k)%DIL)))
      zero_bvs_D_tg(k)%DRL  = 0d0
   ENDDO
   
   
   ! (2)
   ! create the Mass matrix
   !
   ALLOCATE( Mass%i      (SIZE(Lns%i))       ); Mass%i       = Lns%i
   ALLOCATE( Mass%i_mumps(SIZE(Lns%i_mumps)) ); Mass%i_mumps = Lns%i_mumps
   ALLOCATE( Mass%j      (SIZE(Lns%j))       ); Mass%j       = Lns%j
   ALLOCATE( Mass%e      (SIZE(Lns%e))       ); Mass%e       = 0d0
   
   CALL qc_0y0_zero_sp_M (mm, jj, 1d0, Mass)
   
   ! (2b)
   ! create the Qf = MassV matrix (non-singular mass only for the velocity)

      CALL zEssM (CMPLX(Mass%e,0d0,KIND=8), Mass%j,  Mass%i,  velCmpnnts*np, &
                                                         Qf%e, Qf%j, Qf%i, Qf%i_mumps)
   
      CALL par_mumps_master (INITIALIZATION, 5, Qf, 0)
      CALL par_mumps_master (SYMBO_FACTOR,   5, Qf, 0)
      write(*,*) 'Volume forcing: starting numerical factorization of Qf ...'
      CALL par_mumps_master (NUMER_FACTOR,   5, Qf, 0)
      write(*,*) '   numerical factorization of Qf: done !'

   
   ! (2C)
   ! create the Q = Mass matrix (singular mass only for the velocity)
   !
   ALLOCATE( Q%i      (SIZE(Lns%i))       ); Q%i       = Lns%i
   ALLOCATE( Q%i_mumps(SIZE(Lns%i_mumps)) ); Q%i_mumps = Lns%i_mumps
   ALLOCATE( Q%j      (SIZE(Lns%j))       ); Q%j       = Lns%j
   ALLOCATE( Q%e      (SIZE(Lns%e))       ); Q%e       = CMPLX(Mass%e,0d0,KIND=8)
   
   ! (2d)
   ! enforce boundary conditions on Mass
   !
   IF     ( TRIM(p_in%HG_vol_in) == 'vol' ) THEN
      CALL Dirichlet_rc_M (np, js_Axis, js_D_tg, 1d0,  Mass)     !!! 0d0 before.  If Mass (and derived matrices)
                                                                 !!! are summed to other matrices with b.c. already
                                                                 !!! applied (and diagE .NE. 0d0 ), this will be ok anyway!
   ELSEIF ( TRIM(p_in%HG_vol_in) == 'in'  ) THEN
      CALL Dirichlet_c_M  (np, js_Axis, js_D_tg,       Mass)
   ENDIF
   
   ALLOCATE( Mass_cmplx%i      (SIZE(Lns%i))       ); Mass_cmplx%i       = Lns%i
   ALLOCATE( Mass_cmplx%i_mumps(SIZE(Lns%i_mumps)) ); Mass_cmplx%i_mumps = Lns%i_mumps
   ALLOCATE( Mass_cmplx%j      (SIZE(Lns%j))       ); Mass_cmplx%j       = Lns%j
   ALLOCATE( Mass_cmplx%e      (SIZE(Lns%e))       ); Mass_cmplx%e       = CMPLX(Mass%e,0d0,KIND=8)

!   DEALLOCATE( Mass%i, Mass%i_mumps, Mass%j, Mass%e )
   
   ! (2e )
   ! Bf matrix
   IF ( TRIM(p_in%HG_vol_in) == 'vol' ) THEN
      CALL zZeroM (Qf, np_L, Bf)                       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! correct it -----> !!!!!!!!!!!!!!!!!!!!! 
      CALL Dirichlet_c_3d_M_rhs (np, js_Axis, js_D_tg, 0d0,  Bf)  !----> CORRECTION needed : set homogeneous Dirichlet b.c.s !!!!!!
   ELSEIF ( TRIM(p_in%HG_vol_in) == 'in') THEN
      CALL zZeroM (BBf, np_L, Bf)
!      write(*,*) 'size(BBf%e) , minval , maxval = ' , size(BBf%e), minval(dble(BBf%e)), maxval(dble(BBf%e))
!      write(*,*) 'size(BBf%i) , minval , maxval = ' , size(BBf%i), minval(     BBf%i ), maxval(     BBf%i )
!      write(*,*) 'size(BBf%j) , minval , maxval = ' , size(BBf%j), minval(     BBf%j ), maxval(     BBf%j )
!      write(*,*) 'size(Bf%e ) , minval , maxval = ' , size(Bf%e ), minval(dble(Bf%e )), maxval(dble(Bf%e ))
!      write(*,*) 'size(Bf%i ) , minval , maxval = ' , size(Bf%i ), minval(     Bf%i  ), maxval(     Bf%i  )
!      write(*,*) 'size(Bf%j ) , minval , maxval = ' , size(Bf%j ), minval(     Bf%j  ), maxval(     Bf%j  )
   ENDIF
   
!   write(*,*) 'BBf%j = ' , BBf%j
!   write(*,*) 'Bf%j  = ' , Bf%j
!   STOP

  
   ! (3)
   ! create the matrices for the DIRECT AND ADJOINT PROBLEM. 
   !
   ! Direct problem
   ALLOCATE( Lns_cmplx%i      (SIZE(Lns%i))       ); Lns_cmplx%i       = Lns%i
   ALLOCATE( Lns_cmplx%i_mumps(SIZE(Lns%i_mumps)) ); Lns_cmplx%i_mumps = Lns%i_mumps
   ALLOCATE( Lns_cmplx%j      (SIZE(Lns%j))       ); Lns_cmplx%j       = Lns%j
   ALLOCATE( Lns_cmplx%e      (SIZE(Lns%e))       ); Lns_cmplx%e       = CMPLX(0d0,0d0,KIND=8)

   Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
   CALL extract(x_vec, u0)
   CALL qc_1y1_sp_gg_3d_M  (mm, jj,                 1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
   CALL qc_oseen2y_sp_3d_M (mm, jj,                  u0,    beta,  Lns_cmplx) ! + linearized terms
   CALL qc_1y0_sp_3d_M     (mm, jj, jj_L,          -1d0,    beta,  Lns_cmplx) ! + pressure gradient (ibp)
   CALL qc_0y1_sp_3d_M     (mm, jj, jj_L,          -1d0,    beta,  Lns_cmplx) ! - velocity divergence
   IF ( TRIM(p_in%HG_vol_in) == 'vol' ) THEN
      CALL Dirichlet_rc_3d_M  (np, js_Axis, js_D_tg,   1d0, Lns_cmplx) ! Dirichlet BCs
   ELSEIF ( TRIM(p_in%HG_vol_in) == 'in' ) THEN
      CALL Dirichlet_c_3d_M   (np, js_Axis, js_D_tg,        Lns_cmplx) ! Dirichlet BCs
   END IF

!   write(*,*) ' check point 1'
!   write(*,*) ' check point 2'
!   write(*,*) ' check point 3'
!   write(*,*) ' check point 4'

!   WRITE(*,*) 'desingularize_tg = ' , DESINGULARIZE_tg $$$$$$$$$$    FALSE   $$$$$$$$$$$$
   IF (DESINGULARIZE_tg) THEN								!$
      ! row										!$
      DO i = Lns_cmplx%i(Nx), Lns_cmplx%i(Nx + 1) - 1					!$
          Lns_cmplx%e(i) = CMPLX(0d0,0d0,KIND=8)					!$
         IF (Lns_cmplx%j(i) == Nx) Lns_cmplx%e(i) = CMPLX(1d0,0d0,KIND=8)		!$
      ENDDO										!$
   ENDIF										!$
!write(*,*) 'maxvals0d ', MAXVAL(DBLE(Lns_cmplx%e)), MAXVAL(AIMAG(Lns_cmplx%e))		!$
!write(*,*) 'minvals0d ', MINVAL(DBLE(Lns_cmplx%e)), MINVAL(AIMAG(Lns_cmplx%e))		!$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   ALLOCATE( Wd%i      (SIZE(Lns%i))       ); Wd%i       = Lns%i
   ALLOCATE( Wd%i_mumps(SIZE(Lns%i_mumps)) ); Wd%i_mumps = Lns%i_mumps
   ALLOCATE( Wd%j      (SIZE(Lns%j))       ); Wd%j       = Lns%j
   ALLOCATE( Wd%e      (SIZE(Lns%e))       ); Wd%e       = CMPLX(0d0,0d0,KIND=8)

   CALL par_mumps_master (INITIALIZATION, 6, Wd, 0)
   CALL par_mumps_master (SYMBO_FACTOR  , 6, Wd, 0)

!   write(*,*) ' check point 5'
!   write(*,*) ' check point 6'
!   write(*,*) ' check point 7'
   
   ! Adjoint problem
!!!   Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
!!!   CALL extract(x_vec, u0)
!!!   CALL qc_1y1_sp_gg_3d_M  (mm, jj,                 1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
!!!   CALL qc_oseen2y_sp_3d_M (mm, jj,                  u0,    beta,  Lns_cmplx) ! + linearized terms
!!!   CALL qc_1y0_sp_3d_M     (mm, jj, jj_L,          -1d0,    beta,  Lns_cmplx) ! + pressure gradient (ibp)
!!!   CALL qc_0y1_sp_3d_M     (mm, jj, jj_L,          -1d0,    beta,  Lns_cmplx) ! - velocity divergence
!!!   CALL Dirichlet_rc_3d_M  (np, js_Axis, js_D_tg,   1d0,           Lns_cmplx) ! Dirichlet BCs
!!!!   WRITE(*,*) 'desingularize_tg = ' , DESINGULARIZE_tg $$$$$$$$$$    FALSE   $$$$$$$$$$$$
!!!   IF (DESINGULARIZE_tg) THEN
!!!   ! column
!!!   WHERE (Lns_cmplx%j == Nx)
!!!         Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
!!!   ENDWHERE
!!!   IF (Lns_cmplx%j(SIZE(Lns_cmplx%j,1)) == Nx) Lns_cmplx%e(SIZE(Lns%j,1)) = CMPLX(1d0,0d0,KIND=8)
!!!   ENDIF
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!! INIT and SYM FACTORIZATION OF Wa : INEFFICIENT IMPLEMENTATION OF THE CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <------------------
!!  It is possible to use Wd and the TRANSPOSE solver of MUMPS, w/o computing Wa, if CONJG terms are properly used.
!!   ALLOCATE( Wa%i      (SIZE(Lns%i))       ); Wa%i       = Lns%i
!!   ALLOCATE( Wa%i_mumps(SIZE(Lns%i_mumps)) ); Wa%i_mumps = Lns%i_mumps
!!   ALLOCATE( Wa%j      (SIZE(Lns%j))       ); Wa%j       = Lns%j
!!   ALLOCATE( Wa%e      (SIZE(Lns%e))       ); Wa%e       = CMPLX(0d0,0d0,KIND=8)
!!  
!!   CALL par_mumps_master (INITIALIZATION, 7, Wa, 0)
!!   CALL par_mumps_master (SYMBO_FACTOR  , 7, Wa, 0)
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <------------------
   
! ------------------------------------------------------------------------------
 
   p_in%tranGrowth_initGuess = 1 ! forcing random guess by ARPACK
   SELECT CASE ( p_in%tranGrowth_initGuess )
	 CASE (1)
	    ! (a)
	    ! random guess created by ARPACK
	    WRITE(*,*) '    initial guess: ARPACK''s random'
	    WRITE(*,*)
	 CASE DEFAULT
	    ! (a)
	    ! random guess created by ARPACK
	    WRITE(*,*) '    initial guess: ARPACK''s random'
	    WRITE(*,*)
   END SELECT


! ==============================================================================
!   REAL   (KIND=8), DIMENSION(:)  , ALLOCATABLE :: sumGains
!   REAL   (KIND=8), DIMENSION(:,:), ALLOCATABLE :: sumBaseFlowSens
!   REAL   (KIND=8), DIMENSION(:,:), ALLOCATABLE :: sumVolForceSens
!   REAL   (KIND=8), DIMENSION(:,:), ALLOCATABLE :: sumWallBlowSuctSens
! Allocate variables for the Stochastic gain computation
   IF ( p_in%Stoc_Id == 'yes' ) THEN
      ALLOCATE( sumGains(SIZE(omega)) )
      ALLOCATE( sumBaseFlowSens    (velcmpnnts*np,SIZE(omega)) )
      ALLOCATE( sumVolForceSens    (velcmpnnts*np,SIZE(omega)) )
      ALLOCATE( sumWallBlowSuctSens(velcmpnnts*np,SIZE(omega)) )
      sumGains = 0.0d0
      sumBaseFlowSens     = 0.0d0
      sumVolForceSens     = 0.0d0
      sumWallBlowSuctSens = 0.0d0
   END IF   
! ==============================================================================
 ! cycle on Omega   
   DO ii = 1, SIZE(omega,1)
        WRITE(*,*)
        WRITE(*,*) '*************************************'
        write(*,*) '**** omega(', ii, ') = ' , omega(ii)
        WRITE(*,*) '*************************************'
        WRITE(*,*)


! ==============================================================================
! Parameters setting for ZNEUPD
! ido :  REVERSE COMMUNICATION FLAG: 
  !        primo passo ido = 0; poi parametro di zneupd per dire se 
  !        c'è convergenza o no
    ido = 0
  ! bmat : matrix for the semi-inner product 
  !        'I' for std eigenproblem
  !        'G' for generalized e.p.
    bmat = 'I'
  ! n : dimension of the e.p.
    IF ( TRIM(p_in%HG_vol_in) == 'vol' ) THEN
       n = velcmpnnts * np           !nX
    ELSEIF ( TRIM(p_in%HG_vol_in) == 'in' ) THEN
       n = velcmpnnts * SIZE(nodes_Inflow_Raw)
    ENDIF
    
! which : 'LM' -> want the NEV eigenvalues of largest magnitude.
!          'SM' -> want the NEV eigenvalues of smallest magnitude.
!          'LR' -> want the NEV eigenvalues of largest real part.
!          'SR' -> want the NEV eigenvalues of smallest real part.
!          'LI' -> want the NEV eigenvalues of largest imaginary part.
!          'SI' -> want the NEV eigenvalues of smallest imaginary part.
    which = 'LM'
  ! nev : Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
    nev = p_in%HG_num_compute
  ! tol :
    tol = 0d0 !0d0 
  ! resid :Complex*16 array of length N.  (INPUT/OUTPUT)
!          On INPUT: 
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
    ALLOCATE(resid(n))
  ! ncv :
    ncv = 2 * nev + 2
  ! v : Complex*16 array N by NCV.  (OUTPUT)
!          Contains the final set of Arnoldi basis vectors. 
    ALLOCATE(v(n,ncv))
  ! ldv : 
    ldv = n
  ! iparam :
    maxit = 500
    iparam(1) = 1
    iparam(3) = maxit
    iparam(7) = 1
  ! ipntr : Integer array of length 14.  (OUTPUT)
!     Pointer to mark the starting locations in the WORKD and WORKL
!     arrays for matrices/vectors used by the Arnoldi iteration.
   ! ipntr 
  ! workd : Complex*16 work array of length 3*N.  (REVERSE COMM.)
!          Distributed array to be used in the basi! Arnoldi iteration
!          for reverse communication.  The user should not use WORKD 
!          as temporary workspace during the iteration !!!!!!!!!!
    ALLOCATE(workd(3*n)) 
  ! lworkl : Integer. (INPUT) LWORKL must be at least 3*NCV**2 + 5*NCV
        !lworkl = ncv*(ncv+8)    ! SYMMETRI!     -> DSAUPD
        !lworkl = 3*ncv**2+6*ncv ! NON SYMMETRI! -> DNAUPD
    lworkl = 3*ncv**2+5*ncv  ! SYMMETRI! -> ZNAUPD 
    ! workl : Complex*16 work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
    ALLOCATE(workl(lworkl))  
  ! rwork : Double precision  work array of length NCV (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
    ALLOCATE(rwork(ncv))
  ! info :
    info = 0
    
! ZNAUPD
    ALLOCATE(d(ncv))
    ALLOCATE(lselect(ncv))
    ALLOCATE(workev(2*ncv))
    sigma = CMPLX(0d0,0d0)
   
   	write(*,*) 'start cycle on Omega'

    ! Q , Qf , Bf already defined ----------------------------------------------  
    
    ! NUMERICAL FACTORIZATION of Wd and Wa
    ! Wd -----------------------------------------------------------------------  
      CALL zAlinB_s (CMPLX(0.0,omega(ii),kind=8), Mass_cmplx%e, &
                     CMPLX(1.0,0.0,KIND=8),   Lns_cmplx%e, Wd%e)
      write(*,*) 'L+i*omega*B : starting numerical factorization'
      CALL par_mumps_master (NUMER_FACTOR, 6, Wd, 0)
      write(*,*) 'L+i*omega*B : completed numerical factorization'
!! NUMERICAL FACTORIZATION of Wa : INEFFICIENT IMPLEMENTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <----------------------      
!!    ! Wa -----------------------------------------------------------------------
!!      CALL zAlinB_s (CMPLX(0.0,-omega(ii),KIND=8), Mass_cmplx%e, &
!!                     CMPLX(1.0,0.0,KIND=8),   CONJG(Lns_cmplx%e), Wa%e)
!!      write(*,*) '(L+i*omega*B)H : starting numerical factorization'
!!      CALL par_mumps_master (NUMER_FACTOR, 7, Wa, 0)
!!      write(*,*) '(L+i*omega*B)H : completed numerical factorization'
!! NUMERICAL FACTORIZATION of Wa : INEFFICIENT IMPLEMENTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <----------------------

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$      
! Actions on f for the REVERSE COMMUNICATION LOOP
  ! y1
  ! y2 = Bf * y1
  ! A * y3 = y2
  ! y4 = Q * y3         ! from ( U , P ) to ( U , 0 )
  ! A^H * y5 = y4
  ! y6 = Bf' * y5
  ! Qf * y7 = y6
      
      IF ( TRIM(p_in%HG_vol_in) == 'vol' ) THEN
         ALLOCATE( y1(velCmpnnts*np) , y6(velCmpnnts*np) , y7(velCmpnnts*np) , STAT=IERR)
      ELSEIF ( TRIM(p_in%HG_vol_in) == 'in' ) THEN
         ALLOCATE( y1(velCmpnnts*SIZE(nodes_Inflow_Raw)) , STAT=aERR)
         ALLOCATE( y6(velCmpnnts*SIZE(nodes_Inflow_Raw)) , STAT=aERR)
         ALLOCATE( y7(velCmpnnts*SIZE(nodes_Inflow_Raw)) , STAT=aERR)
      ENDIF      
      ALLOCATE( y2(nX) , y3(nX) , y4(nX) , y5(nX)  )

      i = 0
      DO
         i = i+1
         WRITE(*,*) 'Harmonic Gain, iteration n. ', i

! ZNAUPD ROUTINE -----------------------------------------------------
         CALL znaupd  ( ido  , bmat , n     , which, &                  !
   	                nev  , tol  , resid , ncv,   &                  !
   	                v    , n    , iparam, ipntr, &                  !
   	                workd, workl, lworkl, rwork, info )             !
!---------------------------------------------------------------------

         IF (ido .eq. -1 .or. ido .eq. 1) THEN
! Action on the vector ---------------------------------------------------

	    y1 = workd(ipntr(1):ipntr(1)+n-1)

!            IF ( TRIM(p_in%HG_vol_in) == 'vol' .AND. TRIM(p_in%HG_control_domain) == 'inlet' ) THEN
!               CALL phiControlDomain (rr,p_in%HG_control_domain_blending,y1)
!            END IF

	    CALL zAtimx(y2, Bf%e, Bf%j, Bf%i, y1)
	    y3 = y2
	    CALL par_mumps_master (DIRECT_SOLUTION, 6, Wd, 0, y3)
	    CALL zAtimx (y4, Q%e, Q%j, Q%i, y3)
!!  FOR THE INEFFICIENT IMPLEMENTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <-----------
!!            y5 = y4
!!            CALL par_mumps_master (TRANSP_SOLUTION, 7, Wa, 0, y5)
!!            CALL zAtimx_T(y6, Bf%e, Bf%j, Bf%i, y5)
!!  FOR THE INEFFICIENT IMPLEMENTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <-----------
            y5 = CONJG(y4)            
	    CALL par_mumps_master (TRANSP_SOLUTION, 6, Wd, 0, y5)
	    y5 = CONJG(y5)
	    CALL zAtimx_T(y6, Bf%e, Bf%j, Bf%i, y5)
            y7 = y6
          
!            IF ( TRIM(p_in%HG_vol_in) == 'vol' .AND. TRIM(p_in%HG_control_domain) == 'inlet' ) THEN
!               CALL phiControlDomain (rr,p_in%HG_control_domain_blending,y7)
!            END IF

            IF     ( TRIM(p_in%HG_vol_in) == 'vol' ) THEN
               CALL par_mumps_master (DIRECT_SOLUTION,  5, Qf   , 0, y7)
            ELSEIF ( TRIM(p_in%HG_vol_in) == 'in'  ) THEN
               CALL par_mumps_master (DIRECT_SOLUTION, 14, Qf_In, 0, y7)
            END IF
	 
	    workd(ipntr(2):ipntr(2)+n-1) = y7
!-------------------------------------------------------------------------
!        %------------------------------------------%
!        | L O O P   B A C K to call ZNAUPD  again. |
!        %------------------------------------------%
	    CYCLE
         ELSE IF ( ido .eq. 99) THEN
	       
	    WRITE(*,*) 'Done Harmonic Gain, exit info: ', info
         EXIT
         END IF 
		 
      END DO ! Reverse Communication Loop
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
	
! Convergence or error? ------------------------------------------------
      IF ( info .lt. 0 ) THEN ! There is an error.
         WRITE(*,*) ' Error with ZNAUPD, info = ', info
         WRITE(*,*) ' Check the documentation of ZNAUPD.'
      ELSE ! We made it.
         WRITE(*,*) 'ok'
! Call to ZNEUPD -------------------------------------------------------         
         rvec = .true.
         CALL zneupd ( rvec , 'A'   , lselect, d     , v,      &
                       n    , sigma , workev , bmat  , n,      &
                       which, nev   , tol    , resid , ncv,    &
                       v    , n     , iparam , ipntr , workd,  &
                       workl, lworkl, rwork  , ierr )
         IF ( ierr .ne. 0 ) THEN
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DSEUPD. |
!           %------------------------------------%
!
            WRITE(*,*) ' Error with ZNEUPD, info = ', ierr
            STOP       ' Check the documentation of DSEUPD.'

         ELSE

            ! Number of converged singular values
            nconv =  iparam(5)

            ALLOCATE(singularValues(nconv))
            ALLOCATE(singularVectors(n,nconv))
            ALLOCATE(singularVectorsP(nX))
            ALLOCATE(optimalResponse(velCmpnnts*np,nconv))
            
            ! Results:
            !    Harmonic Gains   = singularValues
            !    Optimal Forcing  = singularVectors
            !    Optimal Response = R * f
            singularValues  = d(1:nconv)
            singularVectors = v(:,1:nconv)
            DO i3 = 1 , nconv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               IF ( TRIM(p_in%HG_vol_in) == 'vol' .AND. TRIM(p_in%HG_control_domain) == 'inlet' ) THEN
                  CALL phiControlDomain (rr,p_in%HG_control_domain_blending,singularVectors(:,i3))
               END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               CALL zAtimx (singularVectorsP, Bf%e, Bf%j, Bf%i, singularVectors(:,i3))                 
               CALL par_mumps_master (DIRECT_SOLUTION, 6, Wd, 0, singularVectorsP)
               optimalResponse(:,i3) = singularVectorsP(1:velcmpnnts*np)
            ENDDO 
            WRITE(*,*) singularValues
            
            IF ( p_in%Stoc_Id == 'yes' ) THEN
               sumGains(ii) = SUM(singularValues)
            END IF

! File HarmonicGains.txt            
            format_string = "(A2,I2.2)"
            WRITE(omegaString,format_string) 'om' , ii
            
!            WRITE(*,*) 'II          = ' , 'A' , II , 'B'
!            WRITE(*,*) 'omegaString = ' , 'A' // omegaString // 'B'
!            WRITE(*,*) 'TRIM(omegaString) = ' , 'A' // TRIM(omegaString) // 'B'
!            WRITE(*,*) 'TRIM(p_in%HG_output_directory) = ' , 'A' // TRIM(p_in%HG_output_directory) // 'B'
!            WRITE(*,*) 'TRIM(p_in%HG_vol_in) = ' , 'A' // TRIM(p_in%HG_vol_in) // 'B'
!            
!            
!            
!            
!            WRITE(*,*)
!            WRITE(*,*) 'A' // TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
!    '/HarmonicGains_' // TRIM(p_in%HG_vol_in) // '_' // TRIM(omegaString) // '.txt' // 'B'
!            WRITE(*,*)
!            
!            
            

! Create folders ex: './HarmonicGain/case/om' , './HarmonicGain/case/om/f' , './HarmonicGain/case/om/u'

  write(*,*) 'mkdir -p '// TRIM(p_in%HG_output_directory) // TRIM(omegaString)
  WRITE(*,*) 'mkdir -p '// TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in)
  WRITE(*,*) 'mkdir -p '// TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // '/f'
  WRITE(*,*) 'mkdir -p '// TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // '/u'
  CALL system( 'mkdir -p --verbose '// TRIM(p_in%HG_output_directory) // TRIM(omegaString) )
  CALL system( 'mkdir -p --verbose '// TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) )
  CALL system( 'mkdir -p --verbose '// TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // '/f' )
  CALL system( 'mkdir -p --verbose '// TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // '/u' )
  
  
  
   OPEN(UNIT=21,FILE = TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
    '/HarmonicGains_' // TRIM(p_in%HG_vol_in) // '_' // TRIM(omegaString) // '.txt' )
     do i2 = 1 , MIN(sizeOmega,nconv)
         write(21,*) singularValues(i2)
     enddo
   close(21)
! Files optForcing_omXX_kXX.vtk in ./HarmonicGain/case/omegaXX/f
   numHGplot = MIN(nconv,p_in%HG_num_compute,p_in%HG_num_plot)
           
   IF ( TRIM(p_in%HG_vol_in) == 'vol' ) THEN
      CALL vtk_plot_eigenvectors (rr, jj,  DBLE(singularVectors(:,1:numHGplot)), &
          TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
          '/f/optForcing_' // TRIM(omegaString) // '_' // 'Real.vtk' )
      CALL vtk_plot_eigenvectors (rr, jj,  AIMAG(singularVectors(:,1:numHGplot)), &
          TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
          '/f/optForcing_' // TRIM(omegaString) // '_' // 'Imag.vtk' )
   ELSE IF ( TRIM(p_in%HG_vol_in) == 'in' ) THEN
      DO i1 = 1 , numHGplot
        ALLOCATE ( uuInlet(3,np) )
        CALL extract_cmplx ( optimalResponse(:,i1) ,  uuInlet )
                 
        write (numEig, "(I2.2)") i1

        CALL write_BVS_INFLOW (1, DBLE (uuInlet), rr, jjs, sides, TRIM(p_in%HG_output_directory) // &
            TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
            '/f/optForcing_' // TRIM(omegaString) // '_eig' // TRIM(numEig) // 'Real1.dat')
        CALL write_BVS_INFLOW (1, AIMAG(uuInlet), rr, jjs, sides, TRIM(p_in%HG_output_directory) // &
            TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
            '/f/optForcing_' // TRIM(omegaString) // '_eig' // TRIM(numEig) // 'Imag1.dat')
        IF (number_of_sides == 10) THEN
           CALL write_BVS_INFLOW (5, DBLE (uuInlet), rr, jjs, sides, TRIM(p_in%HG_output_directory) // &
              TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
              '/f/optForcing_' // TRIM(omegaString) // '_eig' // TRIM(numEig) // 'Real5.dat')
           CALL write_BVS_INFLOW (5, AIMAG(uuInlet), rr, jjs, sides, TRIM(p_in%HG_output_directory) // &
              TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
              '/f/optForcing_' // TRIM(omegaString) // '_eig' // TRIM(numEig) // 'Imag5.dat')
        END IF
        DEALLOCATE ( uuInlet )
      END DO
   ENDIF

! Files optResponse_omXX_kXX.vtk in ./HarmonicGain/case/omegaXX/u
  CALL vtk_plot_eigenvectors (rr, jj,  DBLE(optimalResponse(:,1:numHGplot)), &
       TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
       '/u/optResponse_' // TRIM(omegaString) // '_' // 'Real.vtk' )
  CALL vtk_plot_eigenvectors (rr, jj,  AIMAG(optimalResponse(:,1:numHGplot)), &
       TRIM(p_in%HG_output_directory) // TRIM(omegaString) //'/' // TRIM(p_in%HG_vol_in) // &
       '/u/optResponse_' // TRIM(omegaString) // '_' // 'Imag.vtk' )
                                           
!!!!!!!!!!!!!!!!!!!!!!!! compute the residual norm, not used for now!
!!!!!!!!!!!!!!!!!!!!!!!         CALL dmout(6, nconv, 3, rd, ncv, -6, &
!!!!!!!!!!!!!!!!!!!!!!!                   'Ritz values (Real, Imag) and realtive residuals')

         
         END IF
      END IF
! ----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
