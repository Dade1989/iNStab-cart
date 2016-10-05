MODULE cm_forcing_connectivity
 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module containing the SUBROUTINEs to connect the contributions from the 
! forcing terms expansion in (q,e) to the terms in (a,e)
!
! Last modified : 2016-03-17
!
! 2016-04-18 : v1.0 : cylinder test ok
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! SUBROUTINEs :
!  influenceTerms  ( mm , qq , a , bb , cmb , numCmb )
!  fromQEtoQA      ( a , m , comb )
!  combinations    ( v , comb )
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 USE cm_type_definitions

 IMPLICIT NONE


 CONTAINS


!----------------------------------------------------------------------

 SUBROUTINE influenceTerms ( mm , qq , a , bb , cmb , numCmb )
 ! Find the combinations of indices (cmb) giving nonzero contribution
 ! to the f_mmpp(mm,qq) forcing term in a,e from the forcing term in q,e
 ! of order (a,bb)
 !
 ! INPUTs  :
 !  ( mm , pp )
 !  (  a , bb )
 ! OUTPUTs :
 !  cmb
 
 IMPLICIT NONE
 
 INTEGER , DIMENSION(:) , INTENT(IN) :: mm , qq
 INTEGER                             :: a
 INTEGER , DIMENSION(:) , INTENT(IN) :: bb
 TYPE(forcingInfluence_type) , DIMENSION(:) , POINTER :: cmb
 INTEGER :: numCmb

INTEGER :: nc , np
INTEGER :: io
INTEGER , DIMENSION(:)   , ALLOCATABLE :: vv
INTEGER , DIMENSION(:)   , ALLOCATABLE :: vvnzRaw , tindRaw
INTEGER , DIMENSION(:)   , ALLOCATABLE :: vvnz , tind
TYPE(forcingCombinations_type) , DIMENSION(:) , ALLOCATABLE :: str
INTEGER , DIMENSION(:)   , ALLOCATABLE :: v
INTEGER , DIMENSION(:,:) , ALLOCATABLE :: matComb

INTEGER , DIMENSION(:)   , ALLOCATABLE :: vin
INTEGER , DIMENSION(:,:) , ALLOCATABLE :: indCmbRaw
INTEGER , DIMENSION(:)   , ALLOCATABLE :: s
TYPE(forcingInfluence_type) , DIMENSION(:) , POINTER ::  cmbRaw

INTEGER :: NNZ

INTEGER :: i1 , i2 , i3

nc = SIZE(mm,1)
np = SIZE(qq,1)
! ! Check ----
! WRITE(*,*) nc , np
! ! ----------

! Check the dimension of qq , bb
IF ( SIZE(qq,1) .NE. SIZE(bb,1) ) THEN
  WRITE(*,*) "error : pp,bb must have the same # of elements"
  STOP
END IF

! Check that no negative input exists
DO i1 = 1 , SIZE(mm,1)
  IF ( mm(i1) .LT. 0 ) THEN
    WRITE(*,*) "error : mm(" , i1 ,"is negative! Program stops"
    STOP
  ENDIF
END DO
IF ( a .LT. 0 ) THEN
  WRITE(*,*) "error : a is negative! Program stops"
  STOP
ENDIF
DO i1 = 1 , np
  IF ( qq(i1) .LT.0 ) THEN
    WRITE(*,*) "error : qq(" , i1 ,"is negative! Program stops"
    STOP
  ENDIF
ENDDO
DO i1 = 1 , np
  IF ( bb(i1) .LT. 0 ) THEN
    WRITE(*,*) "error : qq(" , i1 ,"is negative! Program stops"
    STOP
  ENDIF
ENDDO

io = 0
DO i1 = 1 , np
  IF ( qq(i1) - bb(i1) .LT. 0 ) THEN
    io = 1;
  ENDIF
ENDDO

IF ( io .EQ. 1 ) THEN
  ! Do nothing on cmb. Add an output to deal with this case
  ! ...
  numCmb = 0
ELSE 
  ! ...
  ALLOCATE(vv(nc+np))
  vv(   1:nc   ) = mm
  vv(nc+1:nc+np) = qq - bb
  
  ALLOCATE(vvnzRaw(nc+np),tindRaw(nc+np)) 
  nnz = 0
  DO i1 = 1 , nc+np
    IF ( vv(i1) .NE. 0 ) THEN   
      nnz = nnz + 1
      vvnzRaw(nnz) = vv(i1)
      tindRaw(nnz) = i1
    ENDIF
  END DO

  IF ( nnz == 0 ) THEN
    ! Do nothing on cmb. Add an output to deal with this case
    ! ...
    numCmb = 0
  ELSE
    ALLOCATE(vvnz(nnz),tind(nnz))
    vvnz = vvnzRaw(1:nnz)
    tind = tindRaw(1:nnz)

! ! Check ----
! WRITE(*,*) "size(tind) = " , size(tind,1)
! DO i1 = 1 , SIZE(TIND)
!   WRITE(*,*) tind(i1)
! ENDDO
! ! ----------

 
    ALLOCATE(str(nnz)) 
    DO i1 = 1 , nnz
      CALL fromQEtoAE(a,vvnz(i1),str(i1)%mat)
      str(i1)%Nt = SIZE(str(i1)%mat,1)
    ENDDO

    ALLOCATE(v(nnz))
    DO i1 = 1 , nnz
      v(i1) = str(i1)%Nt
    ENDDO
 
    CALL combinations( v , matComb )
   
! ! ----------------
! WRITE(*,*) nc , np
! ! ----------------
    ALLOCATE(vin(nc+np))
    ALLOCATE(indCmbRaw(SIZE(matComb,1),SIZE(matComb,2)))
    ALLOCATE(cmbRaw(SIZE(matComb,1))) 
    ALLOCATE(s(a))    
    numCmb = 0
    DO i1 = 1 , SIZE(matComb,1)
      s = 0
      DO i2 = 1 , nnz 
        s = s + str(i2)%mat(matComb(i1,i2),:)
      ENDDO
      IF ( PRODUCT(s) .EQ. 0 ) THEN
         ! combination not accepted
      ELSE
        numCmb = numCmb + 1
        indCmbRaw(numCmb,:) = matComb(i1,:)
        
        ALLOCATE(cmbRaw(numCmb)%el (a) )
        DO i2 = 1 , a
          ALLOCATE(cmbRaw(numCmb)%el(i2)%mind (nc+np) ) 
          vin = 0
! ! Check ----
! WRITE(*,*) "vin(tind(i3))  tind(i3)         i1        i2       i3      nnz "
! ! ----------
          DO i3 = 1 , nnz
            vin(tind(i3)) = str(i3)%mat(indCmbRaw(numCmb,i3),i2)
! ! Check ----
! WRITE(*,*) VIN(TIND(I3)) , tind(i3) , i1 , i2 , i3 , nnz
! ! ----------
          ENDDO
          cmbRaw(numCmb)%el(i2)%mind = vin
        ENDDO
      ENDIF
    ENDDO
    ALLOCATE(cmb(numCmb))
    cmb = cmbRaw(1:numCmb)

  ENDIF
  
  
ENDIF
 
 
 END SUBROUTINE influenceTerms
 
! !----------------------------------------------------------------------

RECURSIVE SUBROUTINE fromQEtoAE ( a , m , comb )
! Find all the possible combinations of indices a indices whose sum is m
! fromQEtoAE ( 2 , 3 , comb )
! comb =   2   0   0
!          1   1   0
!          1   0   1
!          0   2   0
!          0   1   1
!          0   0   2  

  IMPLICIT NONE
  
  INTEGER , INTENT(IN) :: a , m
  INTEGER , DIMENSION(:,:) , POINTER :: comb
  
  INTEGER , DIMENSION(:,:) , POINTER :: comb1 , combRaw
  INTEGER :: dim1 , dim2
  
  INTEGER :: vSize
  INTEGER :: i1
  
  vSize = 1000
  
!  IF ( ALLOCATED(comb) ) THEN
!    DEALLOCATE(comb) 
!    WRITE(*,*) 'allocated'
!  ENDIF
  
  IF ( a .GT. 1 ) THEN
    
    ALLOCATE(combRaw(vSize,a)) !,combOld(vSize,a))
    combRaw = 0                !, combOld = 0
    
    DO i1 = 0 , m
      CALL fromQEtoAE(a-1,m-i1,comb1)
      IF ( i1 == 0) THEN
        combRaw(1:SIZE(comb1,1),1) = i1
        combRaw(1:SIZE(comb1,1),2:SIZE(comb1,2)+1) = comb1
        dim1    = SIZE(comb1,1)
      ELSE
        combRaw(SIZE(comb1,1)+1:SIZE(comb1,1)+dim1,:) = combRaw(1:dim1,:)
        combRaw(1:SIZE(comb1,1),1) = i1
        combRaw(1:SIZE(comb1,1),2:SIZE(comb1,2)+1)    = comb1
        dim1 = dim1 + SIZE(comb1,1)
      ENDIF
    END DO
  
    ! Extract the actual comb
    ALLOCATE(comb(dim1,a))
    comb = combRaw(1:dim1,:)
       
  ELSE
    ALLOCATE(comb(1,1))
    comb(1,1) = m  
  ENDIF
  
  
END SUBROUTINE fromQEtoAE

!----------------------------------------------------------------------

SUBROUTINE combinations ( v , comb )
! Find all the possible combinations of INTEGER  ...
! Example
! v    = (/ 3 , 2 /)
! comb = 1  1 
!        1  2
!        2  1
!        2  2
!        3  1
!        3  2
  
  IMPLICIT NONE
  
  INTEGER , DIMENSION(:)   , INTENT(IN)  :: v
  INTEGER , DIMENSION(:,:) , ALLOCATABLE :: comb
  
  INTEGER :: n1 , n2
  INTEGER :: repeatedTerms , numberOfRepetitions , termsPerRepetition
  INTEGER :: ind1 , ind2
  
  INTEGER :: i1 , i2 , i3
  
  n2 = SIZE(v)
  n1 = PRODUCT(v)
  
  ALLOCATE(comb(n1,n2))
  
  repeatedTerms       = 1
  numberOfRepetitions = n1
  
  DO i1 = n2 , 1 , -1
    
    termsPerRepetition  = PRODUCT(v(i1:n2))
    numberOfRepetitions = numberOfRepetitions / v(i1)
    
    DO i2 = 1 , v(i1)
      
      DO i3 = 1 , numberOfRepetitions
        
        ind1 = 1 + (i3-1)*termsPerRepetition + (i2-1)*repeatedTerms
        ind2 =     (i3-1)*termsPerRepetition +  i2   *repeatedTerms
        comb(ind1:ind2,i1) = i2
        
      END DO
      
    END DO
    
    repeatedTerms      = repeatedTerms * v(i1)
    
  END DO
  
END SUBROUTINE combinations

!----------------------------------------------------------------------

 
END MODULE cm_forcing_connectivity
