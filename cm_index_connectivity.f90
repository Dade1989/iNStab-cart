MODULE cm_index_connectivity

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module containing the SUBROUTINEs to organise
! . the terms of the expansion
! . the terms of the expansion involved in the forcing terms
!
!
! last modified 2016-03-16
!
! 2016-04-18 v1.0 : cylinder test ok
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! INTERFACEs  :
!  fromMultiIndexToI
!  fromItoMultiIndex
!  
! SUBROUTINEs :
!  fromMultiIndexStrToI  ( mStr , mm , pp , np , ord , iord )
!  fromMultiIndexMatToI  ( mMat , mm , pp , np ,         j1 )
!
!  fromItoMultiIndexStr  ( mStr , ord , iord , np , mm , pp )
!  fromItoMultiIndexMat  ( mMat ,          j , np , mm , pp )
!  
!  fromMIstrToMImat      ( mStr , mMat )
!  fromMImatToMIstr      ( mMat , mStr )
!  
!  allTermsIndices       ( mmax , nc , np , allTerms )
!  termIndices           (    m , nc , np , pSorted , nterms )
!  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

USE cm_type_definitions
USE cm_sort_module

IMPLICIT NONE


INTERFACE fromItoMultiIndex

  MODULE PROCEDURE fromItoMultiIndexStr , &
                   fromItoMultiIndexMat

END INTERFACE fromItoMultiIndex


INTERFACE fromMultiIndexToI

  MODULE PROCEDURE fromMultiIndexStrToI , &
                   fromMultiIndexMatToI

END INTERFACE fromMultiIndexToI



CONTAINS


! -------------------------------------------------------------------------

SUBROUTINE fromMultiIndexStrToI ( mStr , mm , pp , np , ord , iord )
! from multi-Index(mm,pp)   ------> ord , iord in mStr
! (mm,pp) = mStr%order(ord)%indices(iord,:)
  IMPLICIT NONE

  TYPE(powerExpansion) , INTENT(IN) :: mStr
  INTEGER , INTENT(IN) :: np
  INTEGER , DIMENSION(:) , INTENT(IN) :: mm , pp
  INTEGER :: ord , iord

  INTEGER , DIMENSION(:) , ALLOCATABLE :: vv
  INTEGER :: n , nc
  INTEGER :: io
  INTEGER :: i1 , i2 , i3

  n = SIZE(mStr%order(1)%indices,2)
  nc = n - np
  
  ALLOCATE(vv(n))
  vv(   1:nc   ) = mm
  vv(nc+1:nc+np) = pp
  
  DO i1 = 1 , mStr%maxOrder
    DO i2 = 1 , mStr%order(i1)%nterms
      io = 0
      DO i3 = 1 , n
        IF ( mStr%order(i1)%indices(i2,i3) .NE. vv(i3) ) THEN
          io = 1
        END IF
      END DO
      IF ( io .EQ. 0 ) THEN
        ord  = i1
        iord = i2 
        RETURN
      END IF 
    END DO
  END DO
  
!   ! If the program gets here, it has not found the (mm,pp)
!   WRITE(*,*)
!   WRITE(*,*) " error : the multi-Index (mm,pp) has been not found! "
!   WRITE(*,*) "   ----> PROGRAM STOPS !                             "
!   WRITE(*,*)
!   STOP
  
END SUBROUTINE fromMultiIndexStrToI

! -------------------------------------------------------------------------

SUBROUTINE fromMultiIndexMatToI ( mMat , mm , pp , np , j )
! from  multi-Index(mm,pp) ------> index j in mMat 
! (mm,pp) = mMat(j,:)

IMPLICIT NONE

  INTEGER , DIMENSION(:,:) , INTENT(IN) :: mMat
  INTEGER , INTENT(IN) :: np
  INTEGER , DIMENSION(:) , INTENT(IN) :: mm , pp
  INTEGER :: j

  INTEGER , DIMENSION(:) , ALLOCATABLE :: vv
  INTEGER :: n1 , n2 , nc
  INTEGER :: io
  INTEGER :: i1 , i2 , i3

  n1 = SIZE(mMat,1)
  n2 = SIZE(mMat,2)
  nc = n2 - np
  ALLOCATE(vv(n2))
  vv(   1:nc   ) = mm
  vv(nc+1:nc+np) = pp
 
  j = 0 
  DO i1 = 1 , n1
    io = 0
    DO i2 = 1 , n2
      IF ( mMat(i1,i2) .NE. vv(i2) ) THEN
        io = 1
      END IF
    END DO
    IF ( io .EQ. 0 ) THEN
      j = i1
      RETURN
    END IF 
  END DO

!   ! If the program gets here, it has not found the (mm,pp)
!   WRITE(*,*)
!   WRITE(*,*) " error : the multi-Index (mm,pp) has been not found! "
!   WRITE(*,*) "   ----> PROGRAM STOPS !                             "
!   WRITE(*,*)
!   STOP

END SUBROUTINE fromMultiIndexMatToI

! -------------------------------------------------------------------------

SUBROUTINE fromItoMultiIndexStr ( mStr , ord , iord , np , mm , pp )
! from ord , iord in mStr   ------> multi-Index(mm,pp)
! (mm,pp) = mStr%order(ord)%indices(iord,:)
  
  IMPLICIT NONE
  
  TYPE(powerExpansion) , INTENT(IN) :: mStr
  INTEGER , INTENT(IN) :: ord , iord
  INTEGER , INTENT(IN) :: np
  INTEGER , DIMENSION(:) , ALLOCATABLE :: mm , pp
! INTEGER , DIMENSION(:) , ALLOCATABLE , INTENT(OUT) :: mm , pp
  
  INTEGER :: n , nc
  
  INTEGER :: i1
  
  n  = SIZE(mStr%order(1)%indices,2)
  nc = n - np

  ! to avoid INTENT(OUT) in Declaration of mm , pp
  IF ( ALLOCATED(mm) ) THEN
    DEALLOCATE(mm)
    ALLOCATE(mm(nc))  
  END IF
  IF ( ALLOCATED(pp) ) THEN
    DEALLOCATE(pp)
    ALLOCATE(pp(np))  
  END IF

  !ALLOCATE(mm(nc),pp(np))
  mm = mStr%order(ord)%indices(iord,   1:nc   )
  pp = mStr%order(ord)%indices(iord,nc+1:nc+np)
  
END SUBROUTINE fromItoMultiIndexStr

! -------------------------------------------------------------------------

SUBROUTINE fromItoMultiIndexMat ( mMat ,  j , np , mm , pp )
! from index j in mMat   ------> multi-Index(mm,pp)
! (mm,pp) = mMat(j,:)

  IMPLICIT NONE

  INTEGER , INTENT(IN) , DIMENSION(:,:) :: mMat
  INTEGER , INTENT(IN) :: j
  INTEGER , INTENT(IN) :: np
  INTEGER , DIMENSION(:) , ALLOCATABLE :: mm , pp
! INTEGER , DIMENSION(:) , ALLOCATABLE , INTENT(OUT) :: mm , pp
  
  INTEGER :: n , nc
  
  n  = SIZE(mMat,2)
  nc = n - np
  
  ! to avoid INTENT(OUT) in the declaration of mm , pp
  IF ( ALLOCATED(mm) ) THEN
    DEALLOCATE(mm)
    ALLOCATE(mm(nc))  
  END IF
  IF ( ALLOCATED(pp) ) THEN
    DEALLOCATE(pp)
    ALLOCATE(pp(np))  
  END IF
  
  !ALLOCATE(mm(nc),pp(np))
  mm = mMat(j,   1:nc   )
  pp = mMat(j,nc+1:nc+np)
  
  
END SUBROUTINE fromItoMultiIndexMat

! -------------------------------------------------------------------------

SUBROUTINE fromMImatToMIstr ( mMat , mStr )
  ! Input  :
  !  mMat  : multiIndex matrix ( m1 , m2 ) = ( nc+np , totTerms)
  ! Output :
  !  mStr  : multiIndex structure

IMPLICIT NONE

INTEGER , DIMENSION(:,:) , INTENT(IN) :: mMat
TYPE(powerExpansion) :: mStr

INTEGER :: m1 , m2
INTEGER , DIMENSION(:) , ALLOCATABLE :: nTerm1ord
INTEGER :: ord , ordOld , maxOrd

INTEGER :: i1 , j

m1 = SIZE(mMat,1)
m2 = SIZE(mMat,2)

maxOrd = SUM(mMat(m1,:))

mStr%totTerms = m1
mStr%maxOrder = maxOrd

ALLOCATE(mStr%order (maxOrd) )

! Find the number of terms for each order
ALLOCATE(nTerm1ord (maxOrd) )
nTerm1ord = 0

DO j = 1 , m1
  i1 = SUM(mMat(j,:)) 
  nTerm1ord(i1) = nTerm1ord(i1) + 1
END DO

! Allocate structures
DO i1 = 1 , maxOrd

  mStr%order(i1)%order  = i1
  mStr%order(i1)%nterms = nTerm1ord(i1)

  ALLOCATE( mStr%order(i1)%indices (nTerm1ord(i1),m2) )
  mStr%order(i1)%indices = 0

END DO


! Build the structure
ord = 0
i1 = 0
DO j = 1 , m1
  
  ordOld = ord
  ord = SUM(mMat(j,:))

  IF  ( ordOld .EQ. ord ) THEN
    i1 = i1 + 1
  ELSE
    i1 = 1
  ENDIF
  
  mStr%order(ord) % indices(i1,:) = mMat(j,:) 
  
END DO



END SUBROUTINE fromMImatToMIstr

! -------------------------------------------------------------------------

SUBROUTINE fromMIstrToMImat ( mStr , mMat )
  ! Input  :
  !  mStr  : multiIndex structure
  ! Output :
  !  mMat  : multiIndex matrix


  IMPLICIT NONE
  
  TYPE(powerExpansion) , INTENT(IN) :: mStr
  INTEGER , DIMENSION(:,:) , ALLOCATABLE :: mMat
  INTEGER :: nterms
  
  INTEGER :: i1 , i2 , indRow
  
  nterms = mStr%totTerms

  ALLOCATE(mMat ( nterms , SIZE(mStr%order(1)%indices,2) ) )
  
  indRow = 0
  DO i1 = 1 , mStr%maxOrder
    
    DO i2 = 1 , mStr%order(i1)%nterms

      indRow = indRow + 1
      mMat(indRow,:) = mStr%order(i1)%indices(i2,:)
    END DO
    
  END DO 
  
END SUBROUTINE fromMIstrToMImat

! -------------------------------------------------------------------------

SUBROUTINE allTermsIndices ( mmax , nc , np , allTerms )
  ! Input  :
  !  mmax  : order of the expansion
  !  nc    : dimension of the critical space
  !  np    : dimension of the space of the parameters
  ! Output :
  !  pAll : array (nterms,np+nc) containing the indices of all the terms
  !         of order <= m (sorted as ...)
  !  nAll : number of terms of order m
  !  indT : indT(k,1:2) indices of the first and last row in pAll for terms
  !         of order k
  
  IMPLICIT NONE
   
  INTEGER , INTENT(IN) :: mmax , nc , np
  TYPE(powerExpansion) :: allTerms
  INTEGER :: totNum
  
  INTEGER :: i1 , i2
  
  allTerms%maxOrder = mmax

  ALLOCATE(allTerms%order(mmax))
  
  totNum = 0
  DO i1 = 1 , mmax
!! Test -------------
!    WRITE(*,*) "I1 = " , I1 
!    ALLOCATE(allTerms%order(i1)%indices (i1,2) )
!    allTerms%order(i1)%nterms = i1
!    DO i2 = 1 , i1
!      allTerms%order(i1)%indices(:,i2) = i2 
!    END DO
!! End Test ---------

    CALL termIndices ( i1 , nc , np , &
                  allTerms%order(i1)%indices , allTerms%order(i1)%nterms)
    allTerms%order(i1)%order = i1
    totNum = totNum + allTerms%order(i1)%nterms

    
  END DO
  
  allTerms%totTerms = totNum
   
END SUBROUTINE allTermsIndices

! -------------------------------------------------------------------------

SUBROUTINE termIndices ( m , nc , np , pSorted , nterms )
 ! Input  :
 !  m  : order of the terms
 !  nc : dimension of the critical space
 !  np : dimension of the space of the parameters
 ! Output :
 !  pSorted : array (nterms,np+nc) containing the indices of all the terms
 !            of order m (sorted as ...)
 !  nterms  : number of terms of order m
 
 IMPLICIT NONE
 
 INTEGER , INTENT(IN) :: m , nc , np
 INTEGER , DIMENSION(:,:) , POINTER :: pSorted
 INTEGER :: nterms
 
 INTEGER , DIMENSION(:,:) , ALLOCATABLE :: mat , matold
 INTEGER :: hSize
 INTEGER :: m1 , m2  ! Changing dimensions of mat,matold (DO LOOP in 1.)
 
 INTEGER , DIMENSION(:,:), ALLOCATABLE :: po , pRaw , p
 INTEGER :: n , id
 
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: sumNp 
 INTEGER , DIMENSION(:,:), ALLOCATABLE :: pExt , pExtD
 INTEGER                               :: basE
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: basV
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: vSort , iSort
 
 INTEGER :: i1 , i2 , i3
 
 ! 1. Preliminary stage 
 ! Allocate matrices mat, matold bigger than required to avoid repeated
 !  allocation and deallocation. If needed, increase the parameters hSize
 hSize = (nc+np)**m

 ALLOCATE(mat(hSize,m))
 mat    = 0
 ALLOCATE(matold(hSize,m-1))
 matold = 0
 
 DO i1 = 1 , m
   IF ( i1 == 1 ) THEN
     m1 = np+nc
     mat(1:m1,1) = (/ ( i1 , i1 = 1,m1 ) /)
   ELSE
     matold(1:m1,1:m2) = mat(1:m1,1:m2)
     DO i2 = 1 , np + nc
       mat((i2-1)*m1+1:i2*m1,1)      = i2
       mat((i2-1)*m1+1:i2*m1,2:m2+1) = matold(1:m1,1:m2) 
     ENDDO
     m1 = m1 * (np+nc)
   ENDIF
   m2 = i1
 ENDDO
! WRITE(*,*) m1 , (nc+np)**m  

 ! 2. Find indices combinations
 ALLOCATE(po(m1,nc+np))
 po = 0
 DO i1 = 1 , m1
   DO i2 = 1 , m
     po(i1,mat(i1,i2)) = po(i1,mat(i1,i2)) + 1
   ENDDO
 ENDDO

 DEALLOCATE(mat,matold)
 
 ALLOCATE(pRaw(SIZE(po,1),SIZE(po,2)))
 pRaw = 0
 pRaw(1,:) = po(1,:)
 n = 1
 DO i1 = 2 , m1
   id = 0
   DO i2 = 1 , n
     DO i3 = 1 , nc+np
       IF ( pRaw(i2,i3).NE. po(i1,i3) ) THEN
         EXIT
       ENDIF
       IF ( i3 .EQ. nc+np ) THEN
         id = 1
       ENDIF
     ENDDO
   ENDDO
   IF ( id == 0 ) THEN
     n = n+1
     pRaw(n,:) = po(i1,:)
   ENDIF
 ENDDO
 
 ALLOCATE(p(n,nc+np))
 p = pRaw(1:n,:)
 
 nterms  = n

 
 ! 3. Sort indices combinations
 ! 3.a Extended matrices
 ALLOCATE(sumNp(nterms))
 DO i1 = 1 , nterms
   sumNp(i1) = SUM(p(i1,nc+1:nc+np))
 ENDDO
 ALLOCATE(pExt(SIZE(p,1),nc+np+1))
 pExt(:,1:nc+np) = p
 pExt(:,nc+np+1) = sumNp

 basE = m
 ALLOCATE(basV(np+nc+1))
 DO i1 = 1 , nc
   basV(i1) = basE**(nc-i1)
 ENDDO
 DO i1 = 1 , np
   basV(nc+i1) = basE**(np+nc-i1)
 ENDDO
 basV(nc+np+1) = basE**(np+nc)
 
 ALLOCATE(pExtD(SIZE(p,1),nc+np+1))
 pExtD = pExt
 pExtD(:,nc+np+1) = MAXVAL(pExtD(:,nc+np+1)) - pExtD(:,nc+np+1)

 ! 3.b Sorting
 ALLOCATE(vSort(nterms))
 DO i1 = 1 , nterms
   vSort(i1) = SUM( pExtD(i1,:) * basV(:) )
 ENDDO
 
 CALL sort('descend',vSort,iSort)
 ALLOCATE(pSorted(nterms,nc+np))
 pSorted = p(iSort,:)
 
 
END SUBROUTINE termIndices

! -------------------------------------------------------------------------


END MODULE cm_index_connectivity
