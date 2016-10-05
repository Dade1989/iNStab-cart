MODULE hg_arpack_module

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module setting the parameters for ARPACK routines ZNAUPD, ZENUPD for
!  the eigenvalue problem in Harmonic Gain computation
!
! 2016-05-18
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SUBROUTINEs :
!  compute_hg     (  niter , nF , Wd , id_Wd , Q , Bf , Qf , id_Qf )
!  set_ARPACK_par ( fflag , par , nF ) 
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  USE global_variables
  USE prep_mesh_p1p2_sp
  USE par_solve_mumps
  USE sparse_matrix_operations
  USE cartesian_boundary_values
  USE vtk_plot
  USE fem_miscellaneous
  
  IMPLICIT NONE


  
  TYPE arpack_par_type
   ! ZNAUPD -----------------------------------------------------
     INTEGER :: ido , ncv , nev , lworkl , info , n , ldv
     CHARACTER(LEN=1) :: bmat
     CHARACTER(LEN=2) :: which
     INTEGER, DIMENSION(11) :: iparam
     INTEGER, DIMENSION(14) :: ipntr
     REAL(KIND=8) :: tol
     COMPLEX(KIND=8) , DIMENSION(:)  , ALLOCATABLE :: resid
     COMPLEX(KIND=8) , DIMENSION(:)  , ALLOCATABLE :: workd , workl 
     COMPLEX(KIND=8) , DIMENSION(:,:), ALLOCATABLE :: v
     REAL(KIND=8)    , DIMENSION(:)  , ALLOCATABLE :: rwork
   ! ZNEUPD parameters ------------------------------------------
     LOGICAL, DIMENSION(:), ALLOCATABLE :: lselect
     COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: d
     COMPLEX(KIND=8) :: sigma
     COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE :: workev
     INTEGER :: ierr
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: rd
     LOGICAL :: rvec   
     INTEGER :: maxit
  
  END TYPE arpack_par_type
  
  
  CONTAINS
  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE     compute_hg (  niter , nF , &
                             Wd , id_Wd , Q , Bf , &
                             Qf , id_Qf ,          &
                             opt_g , opt_f , opt_u   )

! -------------------------------------------------------------------
! call to ARPACK routines to solve the SINGULAR VALUE problem
! ( u = R f, R is the operator 
! - The singular value problem             :  G^2 = (Rf,Rf)/(f,f)
!   is first tranformed into an EV problem :  R'*R f = G^2 f
! - R'*R f_k = G_k^2 f_k:   G_k         : harmonic gain
!                           f_k         : optimal forcing
!                           u_k = R f_k : optimal response
! - call to ZNAUPD () : cycle until convergence
! - call to ZNEUPD () : reverse loop comunication to get the results
!                       if convergence is reached
! -------------------------------------------------------------------

INTEGER                        , INTENT(IN) :: niter
INTEGER                        , INTENT(IN) :: nF   ! dof of the  
TYPE(CSR_MUMPS_Complex_Matrix) , INTENT(IN) :: Wd , Q , Bf, Qf
INTEGER                        , INTENT(IN) :: id_Wd , id_Qf
COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: opt_g ! HG
COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: opt_f ! opt forcing
COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: opt_u ! opt response 

TYPE(arpack_par_type) :: par

COMPLEX(KIND=8) , DIMENSION(:) , ALLOCATABLE :: y1,y2,y3,y4,y5,y6,y7 

INTEGER :: nconv
COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: opt_fP
COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: opt_f_ex 

CHARACTER(LEN=24) :: format_string
CHARACTER(LEN=24) :: omegaString , eigStr

INTEGER :: IERR
INTEGER :: i1 , i2 , i3 


! set arpack parameters
CALL set_ARPACK_par ( p_in%hg_forcing_flag , par , nF )

! Allocate vctr used in znaupd loop
ALLOCATE( y2(nX) , y3(nX) , y4(nX) , y5(nX)  )

ALLOCATE( y1(nF) , &    ! case dependent : vol
          y6(nF) , &
          y7(nF) , STAT=IERR )
!  ELSEIF ( TRIM(p_in%HG_vol_in) == 'in' ) THEN
!     ALLOCATE( y1(velCmpnnts*SIZE(nodes_Inflow_Raw)) , STAT=aERR)
!     ALLOCATE( y6(velCmpnnts*SIZE(nodes_Inflow_Raw)) , STAT=aERR)
!     ALLOCATE( y7(velCmpnnts*SIZE(nodes_Inflow_Raw)) , STAT=aERR)
!  ENDIF  

! ZNAUPD loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
i1 = 0
DO 
  i1 = i1 + 1
  WRITE(*,*) 'Harmonic Gain, iteration n. ', i1

  WRITE(*,*)  "par%info = " , par%info
  IF ( par%info .NE. 0 ) THEN ; STOP ; ENDIF

! ZNAUPD ROUTINE -----------------------------------------------------
  CALL znaupd  ( par%ido  , par%bmat , par%n     , par%which, &
                 par%nev  , par%tol  , par%resid , par%ncv,   &
                 par%v    , par%n    , par%iparam, par%ipntr, &
                 par%workd, par%workl, par%lworkl, par%rwork, &
                 par%info )

  IF (par%ido .eq. -1 .or. par%ido .eq. 1) THEN

! Action on the vector -----------------------------------------------
    y1 = par%workd(par%ipntr(1):par%ipntr(1)+par%n-1)
    CALL zAtimx(y2, Bf%e, Bf%j, Bf%i, y1)
    y3 = y2
    CALL par_mumps_master (DIRECT_SOLUTION, id_Wd, Wd, 0, y3)
    CALL zAtimx (y4, Q%e, Q%j, Q%i, y3)
    y5 = CONJG(y4)            
    CALL par_mumps_master (TRANSP_SOLUTION, id_Wd, Wd, 0, y5)
    y5 = CONJG(y5)
    CALL zAtimx_T(y6, Bf%e, Bf%j, Bf%i, y5)
    y7 = y6
    CALL par_mumps_master (DIRECT_SOLUTION, id_Qf, Qf, 0, y7)
 
    par%workd(par%ipntr(2):par%ipntr(2)+par%n-1) = y7

    CYCLE
!   Convergence on the required e.v. is not reached yet:
!   --->   L O O P   B A C K to call ZNAUPD  again.

  ELSE IF ( par%ido .eq. 99) THEN
    WRITE(*,*) 'Done Harmonic Gain, exit info: ', par%info
    EXIT
!   Convergence
  END IF 
END DO 

! CALL TO ZNEUPD : Reverse Communication Loop ++++++++++++++++++++++++
! Convergence or error? ----------------------------------------------
IF ( par%info .lt. 0 ) THEN ! There is an error.
  WRITE(*,*) ' Error with ZNAUPD, info = ', par%info
  WRITE(*,*) ' Check the documentation of ZNAUPD.'
ELSE ! We made it.
  WRITE(*,*) 'ok'
! Call to ZNEUPD -----------------------------------------------------         
  par%rvec = .true.
  CALL zneupd ( par%rvec , 'A'        , par%lselect, par%d      , & 
                par%v    , par%n      , par%sigma  , par%workev , &
                par%bmat , par%n      , par%which  , par%nev    , &
                par%tol  , par%resid  , par%ncv    , par%v      , &
                par%n    , par%iparam , par%ipntr  , par%workd  , &
                par%workl, par%lworkl , par%rwork  , par%ierr )
  IF ( par%ierr .ne. 0 ) THEN
  !  Error condition: Check the documentation of ZNEUPD.
    WRITE(*,*) ' Error with ZNEUPD, info = ', ierr
    STOP       ' Check the documentation of DSEUPD.'
  ELSE
 
    ! Number of converged singular values
    nconv =  par%iparam(5)
    WRITE(*,*) " n.conv = " , nconv

    IF ( ALLOCATED(opt_g) ) DEALLOCATE(opt_g)
    IF ( ALLOCATED(opt_f) ) DEALLOCATE(opt_f)
    IF ( ALLOCATED(opt_fP)) DEALLOCATE(opt_fP)
    IF ( ALLOCATED(opt_u) ) DEALLOCATE(opt_u)
    ALLOCATE(opt_g (nconv))
    ALLOCATE(opt_f (par%n,nconv))
    ALLOCATE(opt_fP(nX))
    ALLOCATE(opt_u (velCmpnnts*np,nconv))
    
    ! Results:
    !    Harmonic Gains   = opt_g
    !    Optimal Forcing  = opt_f
    !    Optimal Response = R * f
    opt_g = par%d(1:nconv)
    opt_f = par%v(:,1:nconv)

    DO i2 = 1 , nconv
       CALL zAtimx (opt_fP, Bf%e, Bf%j, Bf%i, opt_f(:,i2)) 
       CALL par_mumps_master (DIRECT_SOLUTION, id_Wd, Wd, 0, opt_fP)
       opt_u(:,i2) = opt_fP(1:velcmpnnts*np)
    ENDDO 
    WRITE(*,*) opt_g
    
!    IF ( p_in%Stoc_Id == 'yes' ) THEN
!       sumGains(ii) = SUM(opt_g)
!    END IF

    ! Plot
    format_string = "(A2,I2.2)"
    WRITE(omegaString,format_string) 'om' , niter

    IF ( p_in%hg_forcing_flag == 'vol' ) THEN
    CALL vtk_plot_eigenvectors (rr, jj, DBLE(opt_f(:,1:5)), &
          TRIM(p_in%hg_output_directory) // TRIM(omegaString) //  &
                                            'fOpt_real.vtk' )
    CALL vtk_plot_eigenvectors (rr, jj,AIMAG(opt_f(:,1:5)), &
          TRIM(p_in%hg_output_directory) // TRIM(omegaString) //  &
                                            'fOpt_imag.vtk' )
    ELSEIF ( p_in%hg_forcing_flag == 'in' ) THEN
      DO i2 = 1 , 5
        CALL extract_z   ( opt_u(:,i2) , opt_f_ex )
        WRITE(eigStr,"(I2.2)") i2
        DO i3 = 1 , SIZE( sides_label_structure%inflow )
          CALL write_BVS ( sides_label_structure%inflow(i3),        & 
                           DBLE(opt_f_ex), rr, jjs, sides,          &
            TRIM(p_in%hg_output_directory) // TRIM(omegaString) //  &
                         'fOpt_real' // "_Eig" // TRIM(eigStr)        )
          CALL write_BVS ( sides_label_structure%inflow(i3),        & 
                           AIMAG(opt_f_ex), rr, jjs, sides,         &
            TRIM(p_in%hg_output_directory) // TRIM(omegaString) //  &
                         'fOpt_imag' // "_Eig" // TRIM(eigStr)        )
        ENDDO
      ENDDO

    ENDIF
 
    CALL vtk_plot_eigenvectors (rr, jj, DBLE(opt_u(:,1:5)), &
          TRIM(p_in%hg_output_directory) // TRIM(omegaString) //  &
                                            'uOpt_real.vtk' )
    CALL vtk_plot_eigenvectors (rr, jj,AIMAG(opt_u(:,1:5)), &
          TRIM(p_in%hg_output_directory) // TRIM(omegaString) //  &
                                            'uOpt_imag.vtk' )

    OPEN(UNIT=21,FILE = TRIM(p_in%HG_output_directory) //   & 
                                     TRIM(omegaString) // 'HG.txt' )
    DO i2 = 1 , nconv
      WRITE(21,*) opt_g(i2)
    ENDDO
    CLOSE(21) 


  END IF
END IF

   
   

END SUBROUTINE compute_hg

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  SUBROUTINE     set_ARPACK_par ( fflag , par , nF )  
    
    IMPLICIT NONE
    
    CHARACTER(*) , INTENT(IN) :: fflag
    TYPE(arpack_par_type)     :: par
    INTEGER , INTENT(IN) , OPTIONAL :: nF

! Check ----
WRITE(*,*) ' velCmpnnts = ' , velCmpnnts
WRITE(*,*) ' np         = ' , np
WRITE(*,*) ' nIn        = ' , nF
! ----------

    par % ido     =  0
    par % bmat    = 'I'
    IF     ( TRIM(fflag) == 'vol' ) THEN
      par % n = velCmpnnts * np
    ELSEIF ( TRIM(fflag) == 'in'  ) THEN
      par % n = nF
    ELSE
      WRITE(*,*) " ERROR in set_ARPACK_par "
      STOP
    ENDIF
    par % which     = 'LM'
    par % nev       = p_in % hg_num_compute
    par % tol       = 0.0d0
    ALLOCATE(par % resid(par % n))
    par % ncv       = 2 * par % nev  + 2
    ALLOCATE(par % v (par%n,par%ncv) )
    par % ldv       = par % n
    par % maxit     = 500
    par % iparam(1) = 1 
    par % iparam(3) = par % maxit 
    par % iparam(7) = 1
    ALLOCATE(par % workd(3 * par % n))
    par % lworkl    = 3*(par%ncv)**2 + 5*par%ncv
    ALLOCATE(par%workl(par%lworkl))
    ALLOCATE(par%rwork(par%ncv))
    par % info      = 0
    ALLOCATE(par%d(par%ncv))
    ALLOCATE(par%lselect(par%ncv))
    ALLOCATE(par%workev (2*par%ncv))
    par % sigma = CMPLX(0.0d0,0.0d0,KIND=8)    ! Shift
  
  END SUBROUTINE set_ARPACK_par

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE hg_arpack_module
