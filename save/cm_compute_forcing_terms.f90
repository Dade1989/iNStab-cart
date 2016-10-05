MODULE cm_compute_forcing_terms

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module containing the SUBROUTINEs to compute the forcing terms
!
!
! 2016-04-18 v1.0: cylinder test ok
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! SUBROUTINEs :
!  computeCm     ( sig , mm , cm )
!  computeFmmpp  ( x0,e0,   QQ,matIglob,      mm,pp,a,bb, f  )
!  computeFmmpp0 (       GG,QQ,matIGlob,nc,np,mm,pp,      fmp)  
!  computeFmmpp2 (       GG,QQ,matIGlob,      mm,pp,      fmp)  
!  
!  !!! check if all the inputs are used !!!
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 USE global_variables
 USE cm_index_connectivity
 USE cm_forcing_connectivity
 USE sparse_matrix_profiles
 USE cm_type_definitions
 USE cm_forcing_input

 IMPLICIT NONE
 
 CONTAINS

! ----------------------------------------------------------------------
 
 SUBROUTINE computeCm( sig , mm , cm )
 
 IMPLICIT NONE
 
 COMPLEX(KIND=8), DIMENSION(:) , INTENT(IN) :: sig
 INTEGER        , DIMENSION(:) , INTENT(IN) :: mm
 COMPLEX(KIND=8) :: cm

 cm = SUM( sig*DBLE(mm) )
 
 END SUBROUTINE computeCm
 
! ----------------------------------------------------------------------

 SUBROUTINE computeFmmpp (x0,e0,js_D,cm_Stiffness_Cmplx,QQ,matIglob, &
                                                           mm,pp,a,bb, f)

 IMPLICIT NONE
 
 REAL(KIND=8)    , DIMENSION(:)   , INTENT(IN) :: x0 , e0
 TYPE(dyn_int_line) , DIMENSION(velCmpnnts) , INTENT(IN) :: js_D
 TYPE(CSR_MUMPS_Complex_Matrix)   , INTENT(IN) :: cm_Stiffness_cmplx
 COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: QQ
 INTEGER         , DIMENSION(:,:) , INTENT(IN) :: matIglob
 INTEGER         , DIMENSION(:)   , INTENT(IN) :: mm , pp 
 INTEGER         , DIMENSION(:,:) , INTENT(IN) :: bb
 INTEGER         , DIMENSION(:)   , INTENT(IN) :: a
 COMPLEX(KIND=8) , DIMENSION(:)   , ALLOCATABLE :: f , f_add
 ! -----------
 TYPE(forcingInfluence_type) , DIMENSION(:) , POINTER :: cmb
 INTEGER :: numCmb
 ! -----------
 INTEGER :: n , nc , np
 ! -----------
 INTEGER         , DIMENSION(:)   , ALLOCATABLE :: mmf , ppf
 INTEGER , DIMENSION(:) , ALLOCATABLE :: indV
 INTEGER :: yq
 ! ----------
 INTEGER :: i1 , i2 , i3
 
 n  = SIZE(x0,1)
 nc = SIZE(mm,1)
 np = SIZE(pp,1)

 IF ( ALLOCATED(f) ) THEN
   DEALLOCATE(f)
 ENDIF
 ALLOCATE(f(n))
 f = 0.0d0

 DO i1 = 1 , SIZE(a,1) 
  WRITE(*,*) ' Forcing Index = ' , i1
! Disp -----
! WRITE(*,*) " a , bb = " , a(i1) , bb(i1,:)
! Disp -----
   CALL influenceTerms ( mm , pp , a(i1) , bb(i1,:) , cmb , numCmb )
! Disp -----
!   WRITE(*,*) "numCmb = " , numCmb
! ----------
   DO i2 = 1 , numCmb 
! Disp -----
!     WRITE(*,*) '   comb : ' , i2 , '/' , numCmb
!     WRITE(*,*) '   SIZE(cmb(i2)%el,1) = ' , SIZE(cmb(i2)%el,1)
! Disp -----
     ALLOCATE(indV(SIZE(cmb(i2)%el,1)))
     indV = 0
     DO i3 = 1 , SIZE(cmb(i2)%el,1)
       mmf = cmb(i2)%el(i3)%mind(1:nc)
       ppf = cmb(i2)%el(i3)%mind(1+nc:nc+np)
       CALL fromMultiIndexToI (matIglob,mmf,ppf,np,yq)
       indV(i3) = yq
! Disp -----
!       WRITE(*,*) '   Term (',i3,') = ' , mmf , ppf
! Disp -----
     ENDDO

     CALL switchForcing (x0,e0,QQ(:,indV),i1,i2,js_D,cm_Stiffness_cmplx,f_add)
     f = f + f_add
     DEALLOCATE(indV)

   ENDDO

   ! Forcing functions with a = 0  and bb = qq
   IF ( (a(i1) .EQ. 0) .AND. ALL(mm .EQ. 0) .AND. ALL(bb(i1,:).EQ.pp) ) THEN
     CALL a0bbForcing ( x0,e0,js_D,cm_Stiffness_cmplx,bb(i1,:), f_add )
     f = f + f_add   
   ENDIF

 ENDDO
 
 END SUBROUTINE computeFmmpp
 
! ----------------------------------------------------------------------

 SUBROUTINE computeFmmpp0 (GG,QQ,matIGlob,nc,np,mm,pp,fmp)
 
 IMPLICIT NONE
 
 COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: GG , QQ
 INTEGER , DIMENSION(:,:) , INTENT(IN) :: matIGlob
 INTEGER , INTENT(IN) :: nc ,np
 INTEGER , DIMENSION(:)   , INTENT(IN) :: mm , pp
 COMPLEX(KIND=8), DIMENSION(:)   , ALLOCATABLE :: fmp
 
 INTEGER , DIMENSION(:)   , ALLOCATABLE :: mm0 , oneS , oneL
 INTEGER :: yg , yq
 INTEGER , DIMENSION(:)   , ALLOCATABLE :: jj , kk 
 
 INTEGER :: i1 , i2

 IF (ALLOCATED(fmp)) THEN
   DEALLOCATE(fmp)
 ENDIF
 ALLOCATE(fmp(SIZE(QQ,1)))
 fmp = 0.0d0
! WRITE(*,*) " SIZE(fmp) = " , SIZE(fmp)

 IF ( SUM(pp) .LE. 0 ) THEN
   RETURN
 ELSE
  
   ALLOCATE(  mm0(SIZE(mm,1)) , oneS(SIZE(pp,1)) )
   ALLOCATE( oneL(SIZE(mm,1)) )
   ALLOCATE(   jj(SIZE(mm,1)) ,   kk(SIZE(pp,1)) )
   
   DO i1 = 1 , np
     mm0  = 0
     oneS = 0
     oneS(i1) = 1
     CALL fromMultiIndexToI (matIGlob,mm0,oneS,np,yg)
     
     IF ( ANY( pp-oneS .GT. 0 ) ) THEN
!     WRITE(*,*) 'io = 1'
       DO i2 = 1 , nc
         oneL = 0
         oneL(i2) = 1
         jj = mm + oneL
         kk = pp - oneS
         CALL fromMultiIndexToI (matIGlob,jj,kk,np,yq)
         IF ( yq .NE. 0) THEN
           IF ( .NOT. ALLOCATED(fmp) ) THEN
             ALLOCATE(fmp(SIZE(QQ,1)))
             fmp = ( mm(i2) + 1) * GG(i2,yg) * QQ(:,yq)
           ELSE
             fmp = fmp + ( mm(i2) + 1) * GG(i2,yg) * QQ(:,yq)
           ENDIF 
         ENDIF
       ENDDO
     ENDIF

   ENDDO
 ENDIF
 
 
 END SUBROUTINE computeFmmpp0

! ----------------------------------------------------------------------

SUBROUTINE computeFmmpp2 (GG,QQ,matIGlob,mm,pp,fmp)
 
 IMPLICIT NONE
 
 COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: GG , QQ
 INTEGER , DIMENSION(:,:) , INTENT(IN) :: matIGlob
 INTEGER , DIMENSION(:)   , INTENT(IN) :: mm , pp
 COMPLEX(KIND=8), DIMENSION(:)   , ALLOCATABLE :: fmp

 INTEGER :: nc ,np
 INTEGER :: ord
 INTEGER , DIMENSION(:,:) , POINTER :: matIndOrd
 INTEGER :: nterm
 
 INTEGER :: yg , yq
 INTEGER , DIMENSION(:)   , ALLOCATABLE :: ii , kk , jj , tt , oneL 
 
 INTEGER ::  i1 , i2 , i3
 
 nc = SIZE(mm,1)
 np = SIZE(pp,1)

 IF (ALLOCATED(fmp)) THEN
   DEALLOCATE(fmp)
 ENDIF
 ALLOCATE(fmp(SIZE(QQ,1)))
 fmp = 0.0d0

 ord = SUM(mm) + SUM(pp)
 IF ( ord .LT. 0 ) THEN
   RETURN
 ELSE
   
   ALLOCATE( oneL(nc) )
   ALLOCATE(   jj(nc) ,   kk(np) )
   ALLOCATE(   ii(nc) ,   tt(np))
   DO i1 = 2 , ord - 1
     CALL  termIndices ( i1 , nc, np , matIndOrd , nterm )
     DO i2 = 1 , nterm
       ii = matIndOrd(i2,   1:nc)
       kk = matIndOrd(i2,nc+1:nc+np)
       CALL fromMultiIndexToI(MatIGlob,ii,kk,np,yq)
       
       DO i3 = 1 , nc
         oneL = 0
         oneL(i3) = 1
         jj = mm - ii + oneL
         tt = pp - kk
         CALL fromMultiIndexToI(MatIGlob,jj,tt,np,yg)
         IF ( yg .NE. 0 ) THEN
           IF ( .NOT. ALLOCATED(fmp) ) THEN
             ALLOCATE(fmp(SIZE(QQ,1)))
             fmp = ii(i3) * GG(i3,yg) * QQ(:,yq)
           ELSE
             fmp = fmp + ii(i3) * GG(i3,yg) * QQ(:,yq)
           ENDIF 
         ENDIF
       ENDDO
     ENDDO
   ENDDO
 ENDIF
 
 

END SUBROUTINE computeFmmpp2

! ----------------------------------------------------------------------

END MODULE cm_compute_forcing_terms
