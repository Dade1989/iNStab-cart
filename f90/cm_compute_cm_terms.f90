MODULE cm_compute_cm_terms

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module containing the SUBROUTINEs to compute
! . the terms of the expansion
!
! 2016-04-18 v1.0 : cylinder test ok
!  - cm_compute_cm_terms : resonant 01 term IF STATEMENT is missing
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! SUBROUTINEs :
!  --- compute10terms ( ... ) TO BE WRITTEN (OR NOT?)
!  compute01terms ( x0,e0,      matIndTot,js_D,L,eL,B,K,PhiC,PsiC,sigC,
!                                                                 a,bb,QQ,GG)
!  computempterms ( x0,e0,order,matIndTot,js_D,L,eL,B,K,PhiC,PsiC,sigC
!                                                                 a,bb,QQ,GG)
!  buildL_cmB     ( L , B , cm  , LcmB )
!  buildExtL_cmB  ( A , B , Phi , Psi , eA ) 
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  USE global_variables
  USE cm_index_connectivity 
  USE cm_compute_forcing_terms
  USE cm_write_complex_arrays
  USE cm_linear
  USE EigenSolve
  USE vtk_plot

  USE stressModule
 
  IMPLICIT NONE
  
  CONTAINS
  
  ! ----------------------------------------------------------------------
  
  ! SUBROUTINE compute10terms ()
  ! 
  ! 
  ! 
  ! END SUBROUTINE compute10terms
  
  ! ----------------------------------------------------------------------
  SUBROUTINE compute01terms ( x0 , e0 , matIndTot, js_D , &
                              L , extL , B , cm_K , PhiC , PsiC , sig , &
                              a , bb , QQ , GG )
  
    IMPLICIT NONE
   
    REAL (KIND=8)   , DIMENSION(:)   , INTENT(IN) :: x0 , e0
    INTEGER         , DIMENSION(:,:) , INTENT(IN) :: matIndTot
    TYPE(dyn_int_line) ,DIMENSION(velCmpnnts) , INTENT(IN) :: js_D
!    REAL   (KIND=8) , DIMENSION(:,:) , INTENT(IN) :: L , B
    TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(IN) :: L , B , extL
    TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(IN) :: cm_K
    COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: PhiC , PsiC
    COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN) :: sig 
    INTEGER , DIMENSION(:)   , INTENT(IN) :: a
    INTEGER , DIMENSION(:,:) , INTENT(IN) :: bb
    COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(INOUT) :: QQ , GG
   
    TYPE(CSR_MUMPS_Complex_Matrix)                 :: Eye_Cmplx
    COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: Phi , Psi
    INTEGER :: ordine
    INTEGER :: nc , np , n
    INTEGER , DIMENSION(:,:) , POINTER :: pSorted
    INTEGER , DIMENSION(:,:) , ALLOCATABLE :: matIndices
    INTEGER :: nterms
    INTEGER , DIMENSION(:) , ALLOCATABLE :: imm , ipp
    ! ------------------
    COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: f01
    COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: B01
    INTEGER , DIMENSION(:) , ALLOCATABLE :: ind0
    ! ------------------
    COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: eigenvalue
    COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: eigenvector
    ! DGESV routine ----
    INTEGER , DIMENSION(:) , ALLOCATABLE :: IPIV
    INTEGER :: INFO
    ! ------------------
    INTEGER :: n0
    COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: AA
    COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: FF
    CHARACTER(LEN=40) :: str
    ! ------------------
    REAL(KIND=8) :: eps
    INTEGER :: term , i1
    
    eps = 1.0e-10
    
    ALLOCATE(Phi(SIZE(PhiC,1),SIZE(PhiC,2))) ; Phi = PhiC
    ALLOCATE(Psi(SIZE(PsiC,1),SIZE(PsiC,2))) ; Psi = PsiC

    n  = SIZE(PhiC,1)
    nc = SIZE(PhiC,2)
    np = SIZE( bb,2)
    ordine = 1
    
    CALL termIndices( ordine , nc , np , pSorted , nterms )
    nterms = nterms - nc
    
    ALLOCATE(matIndices(nterms,SIZE(pSorted,2)))
    matIndices = pSorted(nc+1:nc+nterms,:)
    
    DO term = 1 , nterms
      
      imm = matIndices(term,1:nc)
      ipp = matIndices(term,nc+1:nc+np)
!     WRITE(*,*) "imm = " , imm 
!     WRITE(*,*) "ipp = " , ipp 
      CALL computeFmmpp(x0,e0,js_D,cm_K,QQ,matIndTot,imm,ipp,a,bb,f01)
!      WRITE(*,*) "term = " , term
      WRITE(*,*) "max(abs(f01)) = " , MAXVAL(ABS(f01))


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      ! Check if the smallest eigval of L is smaller than eps
!      WRITE(*,*) "SIZE(L%i) = " , SIZE(L%i)
!      CALL buildEyeMat ( L , Eye_Cmplx )
!      CALL eigensComplexShiftInvert(1                        , &
!                                    p_in%eigen_maxit         , & 
!                                    p_in%eigen_tol           , &
!                                    CMPLX(0.0d0,0.0d0,KIND=8), &
!                                    L         , & ! Jacobian_Cmplx)
!                                    Eye_Cmplx , & ! Eye_Cmplx)
!                                    1   , &
!                                    eigenvalue , &
!                                    eigenvector )
!! Check ----
!WRITE(*,*) "eigenvalue = " , eigenvalue
!
!! ----------
!
!      IF ( MINVAL(ABS(eigenvalue)) .GT. eps ) THEN   ! non-resonant term
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF ( MINVAL(ABS(Sig)) .GT. eps ) THEN   ! non-resonant term
        ! g(:,nc+term) = 0 ;  already set by initialization to ZERO
        WRITE(str,'(A,I1,A,I1,A)') '(A,',nc,'I2,A,',np,'I2,A)'
        WRITE(*,str) "Term (",imm,' |',ipp," ) is non-resonant"
        ! Solve a linear system with matrix L
        CALL par_mumps_master ( NUMER_FACTOR    , 6 , L , 0 )
        
        CALL par_mumps_master ( DIRECT_SOLUTION , 6 , L , 0 , f01 )
        QQ(:,nc+term) = -f01

      ELSE                                 ! resonant term
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! TO BE WRITTEN ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        WRITE(str,'(A,I1,A,I1,A)') '(A,',nc,'I2,A,',np,'I2,A)'
!        WRITE(*,str) "Term (",imm,' |',ipp,") is resonant"
!        n0 = 0
!        ALLOCATE(ind0(nc))
!        DO i1 = 1 , nc
!          IF ( ABS(eigenvalue(i1)) .LT. eps ) THEN
!            n0 = n0 + 1
!            Phi(:,n0) = Phi(:,i1)
!            Psi(:,n0) = Psi(:,i1)
!            ind0(n0) = i1
!          ENDIF
!        ENDDO
!
!        ! Extended Matrix CALL extendMatrix ( L , B , Phi(:,1:n0) , Psi(:,1:n0))
!        ALLOCATE(AA(n+n0,n+n0),FF(n+n0),IPIV(n+n0))
!        AA(1:n,1:n) = L
!        AA(1:n,n+1:n+n0) = -MATMUL(B,Phi(:,1:n0))
!        AA(n+1:n+n0,1:n) =  MATMUL( CONJG(TRANSPOSE(Psi(:,1:n0))),B) 
!        AA(n+1:n+n0,n+1:n+n0) = 0.0d0
!        FF(1:n) = -f01
!        FF(n+1:n+n0) = 0.0d0
!        ALLOCATE(B01(n+n0,1))
!        B01(:,1) = FF
!        CALL CGESV(n+n0,1,AA,n+n0,IPIV,B01,n+n0,INFO)
!        
!        QQ(:,nc+term) = B01(1:n,1)
!        GG(ind0(1:n0),nc+term) = B01(n+1:n+n0,1)
!        DEALLOCATE(AA,FF,B01,IPIV,ind0)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ENDIF

!     CALL writeComplex ( QQ(:,nc+term) , 'QQ' , '(F8.5,A,F8.5,A)' , 8 )
      WRITE(*,*) " nc   = " , nc
      WRITE(*,*) " term = " , term
      CALL writeComplex ( GG(:,nc+term) , 'GG' , '(F8.5,A,F8.5,A)' , 8 )
      WRITE(*,*) "------------------------------------"

      
      CALL vtk_plot_eigenvectors (rr, jj,  &
             DBLE(QQ(:,1:nc+term)), &
             trim(p_in%plot_directory)// &
             'qq001.vtk')

    ENDDO
  
  END SUBROUTINE compute01terms
  
! ----------------------------------------------------------------------

  SUBROUTINE computempterms ( x0 , e0 , order , matIndTot, js_D , &
                              L , extL , B , cm_K , PhiC , PsiC , sig , &
                              a , bb , QQ , GG )
    IMPLICIT NONE
   
    REAL   (KIND=8) , DIMENSION(:)   , INTENT(IN)       :: x0 , e0
    INTEGER                                             :: order
    INTEGER         , DIMENSION(:,:) , INTENT(IN)       :: matIndTot 
    TYPE(dyn_int_line) ,DIMENSION(velCmpnnts) , INTENT(IN) :: js_D
    TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(IN)       :: L , B , extL
    TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(IN)       :: cm_K
    COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN)       :: PhiC , PsiC
    COMPLEX(KIND=8) , DIMENSION(:)   , INTENT(IN)       :: sig 
    INTEGER , DIMENSION(:)   , INTENT(IN) :: a
    INTEGER , DIMENSION(:,:) , INTENT(IN) :: bb
    COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(INOUT) :: QQ , GG

    TYPE(CSR_MUMPS_Complex_Matrix) :: LcmB , eLcmB   
    COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: Phi , Psi
    INTEGER :: nc , np , n
    INTEGER :: y
    INTEGER , DIMENSION(:,:) , POINTER :: pSorted
    INTEGER , DIMENSION(:,:) , ALLOCATABLE :: matIndices
    INTEGER :: nterms
    INTEGER , DIMENSION(:) , ALLOCATABLE :: imm , ipp
    COMPLEX(KIND=8) :: cm
    ! ------------------
    COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: Bf
    COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: fmp , hmp
    COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: eHmp
    COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: h1 , h2 , h3
    COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: B01
    INTEGER , DIMENSION(:) , ALLOCATABLE :: ind0
    ! DGESV routine ----
    INTEGER , DIMENSION(:) , ALLOCATABLE :: IPIV
    INTEGER :: INFO
    ! ------------------
    INTEGER :: n0
    COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: AA
    COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: FF
    CHARACTER(LEN=3)  :: comb_str
    CHARACTER(LEN=40) :: str
    ! ------------------
    REAL(KIND=8) :: eps
    INTEGER :: term , i1 , i2
    
    eps = 1.0e-5
! ! Check ----
! DO i1 = 1 , SIZE(matIndTot,1)
!   WRITE(*,*) matIndTot(i1,:)
! ENDDO
! ! ----------
    
    ALLOCATE(Phi(SIZE(PhiC,1),SIZE(PhiC,2)))
    ALLOCATE(Psi(SIZE(PsiC,1),SIZE(PsiC,2)))

    ALLOCATE(Bf(SIZE(PhiC,1)))

    ALLOCATE(LcmB % i       (SIZE(L%i      ))) ; LcmB%i       = L%i
    ALLOCATE(LcmB % i_mumps (SIZE(L%i_mumps))) ; LcmB%i_mumps = L%i_mumps
    ALLOCATE(LcmB % j       (SIZE(L%j      ))) ; LcmB%j       = L%j
    ALLOCATE(LcmB % e       (SIZE(L%e      ))) ; LcmB%e       = 0.0d0
 
    ALLOCATE(eLcmB % i      (SIZE(extL%i      ))); eLcmB%i       = extL%i
    ALLOCATE(eLcmB % i_mumps(SIZE(extL%i_mumps))); eLcmB%i_mumps = extL%i_mumps
    ALLOCATE(eLcmB % j      (SIZE(extL%j      ))); eLcmB%j       = extL%j
    ALLOCATE(eLcmB % e      (SIZE(extL%e      ))); eLcmB%e       = 0.0d0

    n  = SIZE(PhiC,1)
    nc = SIZE(PhiC,2)
    np = SIZE( bb,2)
    ALLOCATE(hmp(n))


    CALL termIndices( order , nc , np , pSorted , nterms )

    ALLOCATE(matIndices(nterms,SIZE(pSorted,2)))
    matIndices = pSorted
    
    DO term = 1 , nterms
      
      hmp = 0.0d0
      Phi = 0.0d0
      Psi = 0.0d0
      
      imm = matIndices(term,1:nc)
      ipp = matIndices(term,nc+1:nc+np)


      CALL fromMultiIndexToI ( matIndTot , imm , ipp , np , y )

       
      WRITE(*,*) "fmmpp  = ... "
      CALL computeFmmpp(x0,e0,js_D,cm_K,QQ,matIndTot,imm,ipp,a,bb,fmp)
      hmp = -fmp

!+++++
      WRITE(*,*) "fmmpp0 = ... "
      CALL computeFmmpp0 (GG,QQ,matIndTot,nc,np,imm,ipp,fmp)
      CALL zAtimx( Bf , B%e , B%j , B%i , fmp )
      hmp = hmp + Bf !MATMUL(B,fmp)
!+++++

      IF ( order .GT. 2 ) THEN 
        CALL computeFmmpp2 (GG,QQ,matIndTot,imm,ipp,fmp)
        CALL zAtimx( Bf , B%e , B%j , B%i , fmp )
        hmp = hmp + Bf !MATMUL(B,fmp)
      ENDIF
      cm = 0.0d0 
!+++++
      CALL computeCm ( sig , imm , cm )
!+++++

      CALL buildL_cmB ( L , B , cm , LcmB )
      
      IF ( MINVAL(ABS(sig-cm)) .GT. eps ) THEN   ! non-resonant term
      ! g(:,nc+term) = 0 ;  already set by initialization to ZERO
        WRITE(str,'(A,I1,A,I1,A)') '(A,',nc,'I2,A,',np,'I2,A)'
        WRITE(*,str,ADVANCE='no') "Term (",imm,' |',ipp," ) is non-resonant ; "
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL writeComplex( cm , 'cm' , '(F8.5,A,F8.5,A)')

      ! Solve a linear system with matrix L
        CALL par_mumps_master ( NUMER_FACTOR , 6 , LcmB , 0 )
        
        CALL par_mumps_master ( DIRECT_SOLUTION , 6 , LcmB , 0 , hmp )
        QQ(:,y) = hmp

! DGSEV +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  !        WRITE(*,*) " cm = " , cm
!          ALLOCATE(IPIV(n),B01(n,1))
!          B01(:,1) = hmp
!          CALL ZGESV(n,1,L-cm*B,n,IPIV,B01,n,INFO)
!  !        WRITE(*,*) " INFO = " , INFO
!  !        IF (INFO.NE.0) THEN ; WRITE(*,*) "error in CGSEV nr" ; STOP ; ENDIF
!          QQ(:,y) = B01(:,1)
!          DEALLOCATE(IPIV,B01)
! DGSEV +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ELSE                                 ! resonant term
        WRITE(str,'(A,I1,A,I1,A)') '(A,',nc,'I2,A,', np ,'I2,A)'
        WRITE(*,str,ADVANCE='no') "Term (",imm,' |',ipp," ) is resonant     ;  "
        CALL writeComplex( cm , 'cm' , '(F8.5,A,F8.5,A)')
!        WRITE(*,*) " cm = " , cm
        n0 = 0
        ALLOCATE(ind0(nc))
        DO i1 = 1 , nc
          IF ( ABS(SIG(i1)-cm) .LT. eps ) THEN
            n0 = n0 + 1
            Phi(:,n0) = PhiC(:,i1)
            Psi(:,n0) = PsiC(:,i1)
            ind0(n0) = i1
          ENDIF
        ENDDO

! Singular matrix ---------------------
        ALLOCATE(h1(n0)) ; h1 = 0.0d0
        ALLOCATE(h2(n))  ; h2 = 0.0d0
        ALLOCATE(h3(n))  ; h3 = 0.0d0
        DO i1 = 1 , n0
          h1(i1) = -SUM( CONJG(Psi(:,i1))*hmp )
        ENDDO
        GG(ind0(1:n0),y) = h1
        h2 = 0.0d0
        DO i1 = 1 , n0
          h2 = h2 + Phi(:,i1)*h1(i1)
        ENDDO
        CALL zAtimx( h3 , B%e , B%j , B%i , h2 )
        hmp = hmp + h3
        CALL par_mumps_master ( NUMER_FACTOR    , 6 , LcmB , 0 )
        CALL par_mumps_master ( DIRECT_SOLUTION , 6 , LcmB , 0 , hmp )
      ! Projection q = ( I - Phi_c * Psi_c' B ) q --------------
        CALL zAtimx( h2 , B%e , B%j , B%i , hmp )
        DO i1 = 1 , n0
          h1(i1) = -SUM( CONJG(Psi(:,i1))*h2 )
        ENDDO
        h2 = 0.0d0
        DO i1 = 1 , n0
          h2 = h2 + Phi(:,i1)*h1(i1)
        ENDDO
        hmp = hmp + h2
      ! --------------------------------------------------------
        QQ(:,y) = hmp


        DEALLOCATE(ind0)
        DEALLOCATE(h1,h2,h3)




! ! Bordered matrix --------------------
!        CALL buildExtL_cmB ( LcmB , B , Phi(:,1:n0) , Psi(:,1:n0) , eLcmB )
!        ALLOCATE(eHmp(n+n0))
!        eHmp(1:n) = hmp
!        eHmp(n+1:n+n0) = 0.0d0
!
!      ! Solve a linear system with matrix L
!        CALL par_mumps_master ( NUMER_FACTOR , 7 , eLcmB , 0 )
!        
!        CALL par_mumps_master ( DIRECT_SOLUTION , 7 , eLcmB , 0 , eHmp )
!        GG(ind0(1:n0,y)) = eHmp(n+1:n+n0)
!        QQ(:,y) = eHmp(1:n)
!        DEALLOCATE(ind0)
!        DEALLOCATE(eHmp)      GG(ind0(1:n0),y) = eHmp(n+1:n+n0)

! DGSEV +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!          ALLOCATE(AA(n+n0,n+n0),FF(n+n0),IPIV(n+n0))
!          AA(1:n,1:n) = L - cm*B
!          AA(1:n,n+1:n+n0) = -MATMUL(B,Phi(:,1:n0))
!          AA(n+1:n+n0,1:n) =  MATMUL( CONJG(TRANSPOSE(Psi(:,1:n0))),B) 
!          AA(n+1:n+n0,n+1:n+n0) = 0.0d0
!          FF(1:n) = hmp
!          FF(n+1:n+n0) = 0.0d0
!          ALLOCATE(B01(n+n0,1))
!          B01(:,1) = FF
!          CALL ZGESV(n+n0,1,AA,n+n0,IPIV,B01,n+n0,INFO)
!  !        IF (INFO.NE.0) THEN ; WRITE(*,*) "error in CGSEV" ; STOP ; ENDIF
!          
!          QQ(:,y) = B01(1:n,1)
!          GG(ind0(1:n0),y) = B01(n+1:n+n0,1)
!          DEALLOCATE(AA,FF,B01,IPIV,ind0)
! DGSEV +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDIF
!        
!  !     CALL writeComplex ( QQ(:,y) , 'QQ' , '(F8.5,A,F8.5,A)' , 8 )
!      CALL writeComplex ( GG(:,y) , 'GG' , '(F8.5,A,F8.5,A)' , 8 )
      WRITE(*,*) "y       = " , y
      WRITE(*,*) "GG(:,y) = " , GG(:,y)
      WRITE(*,*) "------------------------------------"
      
      WRITE(comb_str,'(3(I1))') imm(1) , imm(2) , ipp
      CALL vtk_plot_eigenvectors (rr, jj,  &
                                  DBLE(QQ(:,y:y)), &
                                  trim(p_in%plot_directory)// &
                                  'qq'// comb_str//'.vtk') 
    ENDDO
  
  END SUBROUTINE computempterms

! ----------------------------------------------------------------------

  SUBROUTINE buildL_cmB ( L , B , cm , LcmB )
  
  IMPLICIT NONE
  
  TYPE(CSR_MUMPS_Complex_Matrix) , INTENT(IN)    :: L , B
  COMPLEX(KIND=8)                , INTENT(IN)    :: cm
  TYPE(CSR_MUMPS_Complex_Matrix) :: LcmB 

  LcmB%e = 0.0d0
  LcmB%e = L%e - cm * B%e
  
  END SUBROUTINE buildL_cmB

! ----------------------------------------------------------------------

  SUBROUTINE buildExtL_cmB ( A , B , Phi , Psi , eA )
  
  IMPLICIT NONE
  
  TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(IN) :: A , B
  COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: Phi , Psi
  
  TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(INOUT) :: eA
  
  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: Bf
  INTEGER :: nn , nPhi , i1
 
  eA%e = 0.0d0 
  nPhi = SIZE(Phi,2)
  nn   = SIZE(Phi,1)
! eA%i , eA%i_mumps , eA%j already defined in start_extended_jacobian
  ALLOCATE( Bf( nn , nPhi  ) )
  DO i1 = 1 , nPHi
    CALL zAtimx( Bf(:,i1) , B%e , B%j , B%i , Phi(:,i1) )
  ENDDO
  DO i1 = 1 , nn
    eA%e(eA%i(i1):eA%i(i1+1)-1-nPhi)   = A%e(A%i(i1):A%i(i1+1)-1)
    eA%e(eA%i(i1+1)-nPhi:eA%i(i1+1)-1) = - Bf(i1,:)
  ENDDO
  DO i1 = 1 , nPhi
    CALL zAtimx_T( Bf(:,i1) , B%e , B%j , B%i , CONJG(Psi(:,i1)) ) 
  ENDDO
  DO i1 = 1 , nPhi 
    eA%e(eA%i(nn+i1):eA%i(nn+i1+1)-1) = Bf(:,i1)
  ENDDO

  END SUBROUTINE buildExtL_cmB

!WRITE(*,*) " Check ------------------"
!WRITE(*,*) " nn   = " , nn
!WRITE(*,*) " nPhi = " , nPhi
!WRITE(*,*) " size(eA%i)       = " , SIZE(eA%i)
!WRITE(*,*) " size(eA%i_mumps) = " , SIZE(eA%i_mumps)
!WRITE(*,*) " size(eA%j)       = " , SIZE(eA%j)
!WRITE(*,*) " size(eA%e)       = " , SIZE(eA%e)
!WRITE(*,*) " END Check --------------"
!
!  WRITE(*,*) nn+nPhi+1
!  WRITE(*,*) eA%i(nn+nPhi+1)
! ----------------------------------------------------------------------

  SUBROUTINE computeCFmpterms ( QQ , side_id , cmCF )

  USE prep_mesh_p1p2_sp
  
  IMPLICIT NONE
  
  COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN)  :: QQ
  INTEGER         , DIMENSION(:)   , INTENT(IN)  :: side_id
  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: cmCF

  COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: qq1
  COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: cf
  COMPLEX(KIND=8) , DIMENSION(:,:) , ALLOCATABLE :: tauW

  INTEGER :: i1  
  
  ALLOCATE(cmCF(2,SIZE(QQ,2)))
  ALLOCATE(qq1(SIZE(QQ,1)))

WRITE(*,*) " Check ---- "
WRITE(*,*) SIZE(qq,2)
WRITE(*,*) " Check ---- "

  OPEN(UNIT=21,FILE='./plots/cmForceCoeff.txt')
  WRITE(21,*) "                    cD                       cL                "
  WRITE(21,*) " ----------- real ------ imag ------- real ------ imag-------- "
  WRITE(21,*) " ------------------------------------------------------------- "
  
  DO i1 = 1 , SIZE(QQ,2)

    qq1 = QQ(:,i1) 
    CALL computeWallStress ( qq1 , rr , side_Id , tauW , cf )
    cmCF(:,i1) = cf
    
    WRITE(21,'(I3,4E17.8)') i1 , DBLE(cf(1)) , AIMAG(cf(1)) , DBLE(cf(2)) , AIMAG(cf(2))

  ENDDO
  
  CLOSE(21)

  
  END SUBROUTINE computeCFmpterms













END MODULE cm_compute_cm_terms
