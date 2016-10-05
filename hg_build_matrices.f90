MODULE hg_build_matrices

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! HARMONIC GAIN
! 
! 2016-05-12
! 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SUBROUTINES : 
!  build_hg_matrices          ( )
!  build_hg_matrices_vol      ( )
!  build_hg_matrices_in       ( ) 
!
!  hgInletVelocityB           ( ) 
!  
!  allocate_mumpsSymbo_direct ( )
!  fill_NumerFact_direct      ( omega , B_e , L_e , Wd , id_Wd )
!
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  USE par_solve_mumps
  USE sparse_matrix_profiles
  USE sparse_matrix_operations
  USE global_variables
  USE prep_mesh_p1p2_sp
  USE qc_sp_M
  USE start_sparse_kit

  USE qb_sp_M
  USE fem_start_sparse
  USE fem_miscellaneous
  USE fem_boundary

  IMPLICIT NONE


  CONTAINS
  
!-----------------------------------------------------------------------
    
  SUBROUTINE build_hg_matrices ( fflag , Lns , x_vec , Mass_c , Lns_c , &
                                 Q , Bf , Qf , id_Qf_mumps , QQ, id_QQ )
    
    IMPLICIT NONE
    
    CHARACTER(*)                , INTENT(IN) :: fflag
    TYPE(CSR_MUMPS_Matrix)      , INTENT(IN) :: Lns
    REAL(KIND=8) , DIMENSION(:) , INTENT(IN) :: x_vec
    TYPE(CSR_MUMPS_Complex_Matrix) :: Mass_c , Lns_c , Q , Bf , Qf , QQ
    INTEGER :: id_Qf_mumps , id_QQ
    
    IF ( TRIM(fflag) == 'vol' ) THEN
      CALL build_hg_matrices_vol ( fflag , Lns , x_vec , Mass_c , Lns_c , &
                                   Q , Bf , Qf , id_Qf_mumps , QQ , id_QQ )
    ELSEIF ( TRIM(fflag) == 'in' ) THEN
      CALL build_hg_matrices_in  ( fflag , Lns , x_vec , Mass_c , Lns_c , &
                                   Q , Bf , Qf , id_Qf_mumps , QQ , id_QQ )
    ELSE
      WRITE(*,*) " ERROR : check the spelling of the variable "
      WRITE(*,*) "         p_in%hg_forcing_flag               "
      STOP
    END IF
  
  END SUBROUTINE build_hg_matrices
    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
  SUBROUTINE build_hg_matrices_vol ( fflag , Lns , x_vec , Mass_c , &
                                     Lns_c , Q , Bf , Qf , id_Qf_mumps , &
                                     QQ , id_QQ )
  
    IMPLICIT NONE
    
    CHARACTER(*)                , INTENT(IN) :: fflag
    TYPE(CSR_MUMPS_Matrix)      , INTENT(IN) :: Lns
    REAL(KIND=8) , DIMENSION(:) , INTENT(IN) :: x_vec
    TYPE(CSR_MUMPS_Complex_Matrix) :: Mass_c , Lns_c , Q , Bf , Qf , QQ
    INTEGER :: id_Qf_mumps , id_QQ

    TYPE(CSR_MUMPS_Matrix) :: Mass , Lns_r

!    TYPE(CSR_MUMPS_Complex_Matrix)      :: Wd
!    INTEGER                             :: id_Wd_mumps

! Volume forcing : +++++++++++++++++++++++++++++++++++++++++++++++++++++
! HOMOGENEOUS boundary conditions must be imposed to Wd * (u,p) = Bf * f
! Wd = i * omega * Mass_c + Lns_c
! SUBROUTINES Dirichlet_rc_M should be used to impose equal to zero both
! the columns and the rows associated to Dirichlet b.c. (nodes js_D); 
! the diagonal element can be set to 1.0d0 (or other values to have a 
! "balanced matrix")

      
  ! Mass matrix ( u+p , real , no b.c. ) : -------------------------------
  !   Mass  ---- 
    ALLOCATE( Mass%i      (SIZE(Lns%i))       ); Mass%i       = Lns%i
    ALLOCATE( Mass%i_mumps(SIZE(Lns%i_mumps)) ); Mass%i_mumps = Lns%i_mumps
    ALLOCATE( Mass%j      (SIZE(Lns%j))       ); Mass%j       = Lns%j
    ALLOCATE( Mass%e      (SIZE(Lns%e))       ); Mass%e       = 0d0
    
    CALL qc_00_zero_sp_M (mm,jj,1.0d0,Mass)
    
  ! Mass matrix ( u+p , cmpl , no b.c. ) : -------------------------------
  !   Q     ---- L2 seminorm (u,p)
    ALLOCATE( Q%i      (SIZE(Lns%i))       ); Q%i       = Lns%i
    ALLOCATE( Q%i_mumps(SIZE(Lns%i_mumps)) ); Q%i_mumps = Lns%i_mumps
    ALLOCATE( Q%j      (SIZE(Lns%j))       ); Q%j       = Lns%j
    ALLOCATE( Q%e      (SIZE(Lns%e)) ); Q%e  = CMPLX(Mass%e,0d0,KIND=8)
    
  ! Mass matrix ( u   , cmpl , no b.c. ) : -------------------------------
  !   Qf    ---- L2 norm (u) or (f)i

    CALL zEssM (CMPLX(Mass%e,0d0,KIND=8), Mass%j,  Mass%i,  velCmpnnts*np, &
                                                Qf%e, Qf%j, Qf%i, Qf%i_mumps)
    id_Qf_mumps = 5
    CALL par_mumps_master (INITIALIZATION, id_Qf_mumps, Qf, 0)
    WRITE(*,*) " Qf :  Starting symfact "
    CALL par_mumps_master (SYMBO_FACTOR,   id_Qf_mumps, Qf, 0)
    WRITE(*,*) " Qf :           symfact done. "
    WRITE(*,*) " Qf :  Starting numfact "
    CALL par_mumps_master (NUMER_FACTOR,   id_Qf_mumps, Qf, 0)
    WRITE(*,*) " Qf :           numfact done. "

  ! Matrix QQ : 'vol' QQ = Qf , 'in' QQ = Qf_vol
    ALLOCATE(QQ%i      (SIZE(Qf%i))       );QQ%i       = Qf%i
    ALLOCATE(QQ%i_mumps(SIZE(Qf%i_mumps)) );QQ%i_mumps = Qf%i_mumps
    ALLOCATE(QQ%j      (SIZE(Qf%j))       );QQ%j       = Qf%j
    ALLOCATE(QQ%e      (SIZE(Qf%e))       );QQ%e       = Qf%e
    id_QQ = id_Qf_mumps


  ! Mass matrix ( u+p , cmpl ,    b.c. ) : -------------------------------
  !   Mass_c
    
    CALL Dirichlet_rc_M (np, js_D, 1d0, Mass) 
     
    ALLOCATE(Mass_c%i      (SIZE(Lns%i))      ); Mass_c%i       = Lns%i
    ALLOCATE(Mass_c%i_mumps(SIZE(Lns%i_mumps))); Mass_c%i_mumps = Lns%i_mumps
    ALLOCATE(Mass_c%j      (SIZE(Lns%j))      ); Mass_c%j       = Lns%j
    ALLOCATE(Mass_c%e (SIZE(Lns%e)) ); Mass_c%e = CMPLX(Mass%e,0d0,KIND=8)
    
  ! From Mass matrix for the RHS ( u+p , cmpl , homo Dir ) : -------------
  !   Bf
    CALL zZeroM (Qf, np_L, Bf)
    CALL Dirichlet_rc_Mz (np, js_D,0.0d0,  Bf)
    !----> set homogeneous Dirichlet b.c.s !!!!!!
     
    DEALLOCATE(Mass%i,Mass%i_mumps,Mass%j,Mass%e)

  ! Jacobian matrix ( u+p , cmpl , b.c ) : -------------------------------    
  ! Lns
!    CALL Dirichlet_rc_M (np, js_D, 0d0,  Lns)
!    ALLOCATE(Lns_c%i      (SIZE(Lns%i))      ); Lns_c%i       = Lns%i
!    ALLOCATE(Lns_c%i_mumps(SIZE(Lns%i_mumps))); Lns_c%i_mumps = Lns%i_mumps
!    ALLOCATE(Lns_c%j      (SIZE(Lns%j))      ); Lns_c%j       = Lns%j
!    ALLOCATE(Lns_c%e  (SIZE(Lns%e)) ); Lns_c%e  = CMPLX(Lns%e ,0d0,KIND=8)
   ALLOCATE( Lns_r%i      (SIZE(Lns%i))       ); Lns_r%i       = Lns%i
   ALLOCATE( Lns_r%i_mumps(SIZE(Lns%i_mumps)) ); Lns_r%i_mumps = Lns%i_mumps
   ALLOCATE( Lns_r%j      (SIZE(Lns%j))       ); Lns_r%j       = Lns%j
   ALLOCATE( Lns_r%e (SIZE(Lns%e)) ); Lns_r%e   = 0.0d0

   ALLOCATE( Lns_c%i      (SIZE(Lns%i))       ); Lns_c%i       = Lns%i
   ALLOCATE( Lns_c%i_mumps(SIZE(Lns%i_mumps)) ); Lns_c%i_mumps = Lns%i_mumps
   ALLOCATE( Lns_c%j      (SIZE(Lns%j))       ); Lns_c%j       = Lns%j
   ALLOCATE( Lns_c%e (SIZE(Lns%e)) ); ! Lns_c%e   = CMPLX(0d0,0d0,KIND=8)

   Lns_r%e = 0d0
   CALL extract(x_vec, u0)
   CALL qc_11_sp_gg_M      (mm, jj,     1d0/Re, Lns_r) ! stifness (GRAD:GRAD)
   CALL qc_oseen2_sp_M     (mm, jj,         u0, Lns_r) ! + linearized terms
   CALL qc_10_sp_M         (mm, jj, jj_L, -1d0, Lns_r) ! + p gradient (ibp)
   CALL qc_01_sp_M         (mm, jj, jj_L, -1d0, Lns_r) ! - velocity div

   CALL Dirichlet_rc_M     (np, js_D,      1d0, Lns_r) ! Dirichlet BCs

   Lns_c % e = CMPLX(Lns_r%e , 0.0d0,KIND=8)


  END SUBROUTINE build_hg_matrices_vol
    
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE build_hg_matrices_in ( fflag , Lns , x_vec , Mass_c , Lns_c , &
                                    Q , Bf , Qf , id_Qf_mumps, QQ, id_QQ  )
  
    IMPLICIT NONE
    
    CHARACTER(*)                , INTENT(IN) :: fflag
    TYPE(CSR_MUMPS_Matrix)      , INTENT(IN) :: Lns
    REAL(KIND=8) , DIMENSION(:) , INTENT(IN) :: x_vec
    TYPE(CSR_MUMPS_Complex_Matrix) :: Mass_c , Lns_c , Q , Bf , Qf , QQ
    INTEGER :: id_Qf_mumps , id_QQ
 
! ---- build Qf
    INTEGER , DIMENSION(:)   , ALLOCATABLE :: mmsIn
    INTEGER , DIMENSION(:,:) , ALLOCATABLE :: jjsIn , jjsInFromOne
    INTEGER , DIMENSION(:)   , ALLOCATABLE :: mms_In
    INTEGER , DIMENSION(:)   , ALLOCATABLE :: nodes_In , nodes_In_Raw
    INTEGER                                :: INFLOW_CENTRE
    TYPE(CSR_MUMPS_Complex_Matrix) :: BBf , qf1
! ----   INTEGER :: id_Qf_mumps

    TYPE(CSR_MUMPS_Matrix) :: Mass , Lns_r

  ! Inlet forcing : ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! INHOMOGENEOUS boundary conditions must be imposed to Wd*(u,p) = Bf * f
  ! the inhomogeneous b.c. is the forcing itself (or something very close)
  ! Wd = i * omega * Mass_c + Lns_c
  ! SUBROUTINES Dirichlet_c_M should be used to impose equal to zero 
  ! the rows associated to Dirichlet b.c. (nodes js_D), except for the 
  ! unitary diagolnal element

      
  ! Mass matrix ( u+p , real , no b.c. ) : -------------------------------
  !   Mass  ---- 
    ALLOCATE( Mass%i      (SIZE(Lns%i))       ); Mass%i       = Lns%i
    ALLOCATE( Mass%i_mumps(SIZE(Lns%i_mumps)) ); Mass%i_mumps = Lns%i_mumps
    ALLOCATE( Mass%j      (SIZE(Lns%j))       ); Mass%j       = Lns%j
    ALLOCATE( Mass%e      (SIZE(Lns%e))       ); Mass%e       = 0d0
    
    CALL qc_00_zero_sp_M (mm,jj,1.0d0,Mass)
    
  ! Mass matrix ( u+p , cmpl , no b.c. ) : -------------------------------
  !   Q     ---- L2 seminorm (u,p)
    ALLOCATE( Q%i      (SIZE(Lns%i))       ); Q%i       = Lns%i
    ALLOCATE( Q%i_mumps(SIZE(Lns%i_mumps)) ); Q%i_mumps = Lns%i_mumps
    ALLOCATE( Q%j      (SIZE(Lns%j))       ); Q%j       = Lns%j
    ALLOCATE( Q%e      (SIZE(Lns%e)) ); Q%e  = CMPLX(Mass%e,0d0,KIND=8)

  ! Mass matrix ( u   , cmpl , no b.c. ) : -------------------------------
  !   QQ    ---- L2 norm (u) or (f)

    CALL zEssM (CMPLX(Mass%e,0d0,KIND=8), Mass%j,  Mass%i,  velCmpnnts*np, &
                                                QQ%e, QQ%j, QQ%i, QQ%i_mumps)
    id_QQ = 5
    CALL par_mumps_master (INITIALIZATION, id_QQ, QQ, 0)
    WRITE(*,*) " Qf :  Starting symfact "
    CALL par_mumps_master (SYMBO_FACTOR,   id_QQ, QQ, 0)
    WRITE(*,*) " Qf :           symfact done. "
    WRITE(*,*) " Qf :  Starting numfact "
    CALL par_mumps_master (NUMER_FACTOR,   id_QQ, QQ, 0)
    WRITE(*,*) " Qf :           numfact done. "
   
  ! Mass matrix for surface elements ( u   , cmpl , no b.c. ) : ----------
  !   Qf 
    CALL boundaryElementsOfSides ( mms , sides , &
                       sides_label_structure%inflow , mmsIn )
    CALL fromJJStoJJSw ( jjs , mmsIn , jjsIn , jjsInFromOne )

    CALL start_matrix_P2SP2S_cmplx     ( JJSIn , qf1 )
    CALL start_matrix_2d_p2SP2S_cmplx  ( JJSInFromOne , qf1 , Qf )

    CALL qs_b_0b0b_s_sp_M  ( mmsIn , jj , jjs , 1.0d0 , Qf)

    id_Qf_mumps = 14
    CALL par_mumps_master (INITIALIZATION, id_Qf_mumps, Qf, 0)
    CALL par_mumps_master (SYMBO_FACTOR,   id_Qf_mumps, Qf, 0)
    write(*,*) 'Inlet forcing: starting numerical factorization of Qf ...'
    CALL par_mumps_master (NUMER_FACTOR,   id_Qf_mumps, Qf, 0)
    write(*,*) '   numerical factorization of Qf: done !'

  ! ~ eye matrix + b.c. at inflow outer nodes
  !    BBf and Bf
  ! nodes_In_Raw (including the extremes) and nodes_In (only "inner" nodes)
    CALL boundaryElementsOfSides(mms,sides,sides_label_structure%inflow,mms_In)
    CALL unionFromMatToVct(jjs(:,mms_In  ),nodes_In_Raw)
      
    CALL inflowNodesOnly(jjs,mms_In,nodes_In_Raw,(/ -1 /),nodes_In,&
                                                            INFLOW_CENTRE)

    CALL hgInletVelocityB ( np , 2 , nodes_In , nodes_In_Raw , BBf )
    CALL zZeroM (BBf, np_L, Bf)

  ! Mass matrix ( u+p , cmpl ,    b.c. ) : -------------------------------
  !   Mass_c
    
    CALL Dirichlet_c_M (np, js_D, Mass) 
     
    ALLOCATE(Mass_c%i      (SIZE(Lns%i))      ); Mass_c%i       = Lns%i
    ALLOCATE(Mass_c%i_mumps(SIZE(Lns%i_mumps))); Mass_c%i_mumps = Lns%i_mumps
    ALLOCATE(Mass_c%j      (SIZE(Lns%j))      ); Mass_c%j       = Lns%j
    ALLOCATE(Mass_c%e (SIZE(Lns%e)) ); Mass_c%e = CMPLX(Mass%e,0d0,KIND=8)
     
    DEALLOCATE(Mass%i,Mass%i_mumps,Mass%j,Mass%e)

  ! Jacobian matrix ( u+p , cmpl , b.c ) : -------------------------------    
  ! Lns
   ALLOCATE( Lns_r%i      (SIZE(Lns%i))       ); Lns_r%i       = Lns%i
   ALLOCATE( Lns_r%i_mumps(SIZE(Lns%i_mumps)) ); Lns_r%i_mumps = Lns%i_mumps
   ALLOCATE( Lns_r%j      (SIZE(Lns%j))       ); Lns_r%j       = Lns%j
   ALLOCATE( Lns_r%e (SIZE(Lns%e)) ); Lns_r%e   = 0.0d0

   ALLOCATE( Lns_c%i      (SIZE(Lns%i))       ); Lns_c%i       = Lns%i
   ALLOCATE( Lns_c%i_mumps(SIZE(Lns%i_mumps)) ); Lns_c%i_mumps = Lns%i_mumps
   ALLOCATE( Lns_c%j      (SIZE(Lns%j))       ); Lns_c%j       = Lns%j
   ALLOCATE( Lns_c%e (SIZE(Lns%e)) ); ! Lns_c%e   = CMPLX(0d0,0d0,KIND=8)

   Lns_r%e = 0d0
   CALL extract(x_vec, u0)
   CALL qc_11_sp_gg_M      (mm, jj,     1d0/Re, Lns_r) ! stifness (GRAD:GRAD)
   CALL qc_oseen2_sp_M     (mm, jj,         u0, Lns_r) ! + linearized terms
   CALL qc_10_sp_M         (mm, jj, jj_L, -1d0, Lns_r) ! + p gradient (ibp)
   CALL qc_01_sp_M         (mm, jj, jj_L, -1d0, Lns_r) ! - velocity div

   CALL Dirichlet_c_M     (np, js_D, Lns_r) ! Dirichlet BCs

   Lns_c % e = CMPLX(Lns_r%e , 0.0d0,KIND=8)
  
  
   
  END SUBROUTINE build_hg_matrices_in


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 SUBROUTINE hgInletVelocityB ( ndof , nblocks , nodesIn , nodesIn_Raw , A )
   
   IMPLICIT NONE 
   
   INTEGER, DIMENSION(:)  , INTENT(INOUT) :: nodesIn
   INTEGER, DIMENSION(:)  , INTENT(INOUT) :: nodesIn_Raw
   INTEGER, INTENT(IN) :: nblocks , ndof
   TYPE(CSR_MUMPS_COMPLEX_Matrix) :: A
   
   INTEGER , DIMENSION(:) , ALLOCATABLE :: NODES
   INTEGER , DIMENSION(:) , ALLOCATABLE :: ORD
   INTEGER :: i1 , i2
   integer :: n_a_d
   
   ! "EXTRA" INFLOW BOUNDARY MODIFICATIONS TO THE BBF IN ORDER TO CONSIDER THE RIGHT B.C.S ON THE AXIS (usually but not
   !     always THE CENTRE OF THE INTAKE MANIFOLD
   ALLOCATE( A%e(nblocks*SIZE(nodesIn_Raw)) ) ;  A%e = CMPLX(1.0d0,0.0d0,KIND=8)
   ALLOCATE( A%i(nblocks*ndof + 1) )    ;
   ALLOCATE( A%j(nblocks*SIZE(nodesIn_Raw)) ) ;
   
   CALL sortDiff  ( nodesIn_Raw , nodes , n_a_d)
   CALL orderSort ( NODES , ord )!!!!!!!! THIS WAY     ORD = 1:size(nodes)
   WRITE(*,*) ORD
   
   DO i1 = 1 , nblocks

      DO i2 = 1 , SIZE(nodes)
         
         IF ( .NOT.(ANY( NODES(I2) == NODESIN ) ) ) THEN
             A%e( (i1-1)*SIZE(nodes) + i2 ) = CMPLX(0.0d0,0.0d0,KIND=8)
         END IF
         
! ! CORRECTION FOR THE AXIALSYMMETRIC CASE: THE CENTRE OF THE AXIS IF THE
!    INFLOW MANIFOLD CONTAINS THE AXIS OF SYMMETRY ...
!        IF (BETA == 0) THEN
!          IF ( (NODES(i2) == INFLOW_CENTRE) .AND. (I1 .NE. 1) ) THEN
!            A%e( (i1-1)*SIZE(nodes) + i2 ) = CMPLX(0.0d0,0.0d0,KIND=8)
!          END IF
!        ELSEIF (BETA == 1) THEN
!          IF ( (NODES(i2) == INFLOW_CENTRE) .AND. (I1 == 1) ) THEN
!            A%e( (i1-1)*SIZE(nodes) + i2 ) = CMPLX(0.0d0,0.0d0,KIND=8)
!          END IF
!        ELSE
!          IF ( NODES(i2) == INFLOW_CENTRE ) THEN
!            A%e( (i1-1)*SIZE(nodes) + i2 ) = CMPLX(0.0d0,0.0d0,KIND=8)
!          END IF
!        END IF
                                                                                            !!! matrix QP-style. Not CSR !!!
         A%j( (i1-1)*SIZE(nodes) + i2 ) = (i1-1)*size(nodes) + ORD(i2)                      !!!!!!!!!!!!!!!!!!!! NO - 1  !!!
         
         A%i( (i1-1)*ndof+1 : (i1-1)*ndof+nodes(1) ) = (i1-1)*SIZE(nodes) + 1               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF ( i2 .LT. SIZE(nodes) ) THEN
            A%i( (i1-1)*ndof + nodes(i2)+1 : (i1-1)*ndof + nodes(i2+1) ) = i2 + (i1-1)*SIZE(nodes) + 1
         ELSE
            A%i( (i1-1)*ndof + nodes(i2)+1 : i1*ndof ) = i2 + (i1-1)*SIZE(nodes) + 1
         END IF
         
      END DO
   
   END DO
   
   A%i( nblocks*ndof + 1 ) = nblocks*SIZE(nodes) + 1
   
   END SUBROUTINE hgInletVelocityB 


! SUBROUTINEs on Wd
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE allocate_mumpsSymbo_direct ( Lns , Wd , id_Wd_mumps )
  
    IMPLICIT NONE

    TYPE(CSR_MUMPS_Matrix) , INTENT(IN) :: Lns 
    TYPE(CSR_MUMPS_Complex_Matrix)      :: Wd
    INTEGER                             :: id_Wd_mumps


    ALLOCATE( Wd%i      (SIZE(Lns%i))       ); Wd%i       = Lns%i
    ALLOCATE( Wd%i_mumps(SIZE(Lns%i_mumps)) ); Wd%i_mumps = Lns%i_mumps
    ALLOCATE( Wd%j      (SIZE(Lns%j))       ); Wd%j       = Lns%j
    ALLOCATE( Wd%e(SIZE(Lns%e))); Wd%e = CMPLX(0d0,0d0,KIND=8)

    id_Wd_mumps = 6
    CALL par_mumps_master (INITIALIZATION, id_Wd_mumps, Wd, 0)
    CALL par_mumps_master (SYMBO_FACTOR  , id_Wd_mumps, Wd, 0)   
 
  END SUBROUTINE allocate_mumpsSymbo_direct

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  SUBROUTINE fill_NumerFact_direct ( omega , B_e , L_e , Wd , id_Wd )
      
    IMPLICIT NONE  
     
    REAL(KIND=8) , INTENT(IN) :: omega 
    COMPLEX(KIND=8) , DIMENSION(:) , INTENT(IN)    :: B_e , L_e
    TYPE(CSR_MUMPS_Complex_Matrix) , INTENT(INOUT) :: Wd
    INTEGER , INTENT(IN) :: id_Wd
    
    
    Wd % e = CMPLX(0.0d0,omega,KIND=8) * B_e + L_e

    WRITE(*,*) ' Wd%i       = ' , SIZE( Wd%i )
    WRITE(*,*) ' Wd%i_mumps = ' , SIZE( Wd%i_mumps )
    WRITE(*,*) ' Wd%j       = ' , SIZE( Wd%j )
    WRITE(*,*) ' Wd%e       = ' , SIZE( Wd%e )
    WRITE(*,*) ' id_wd      = ' , id_Wd
   
    CALL par_mumps_master ( NUMER_FACTOR , id_Wd , Wd , 0 )
     
  END SUBROUTINE fill_NumerFact_direct
  



    
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE hg_build_matrices
