MODULE fem_miscellaneous

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module containing miscellaneous SUBROUTINEs for FEM
!
! Collecting some old routines by Quartapelle and adding other routines
! 2016-05-17
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! INTERFACEs   :
!  fromQPtoCSR
!  fromCSRroQP
!  assembleCSRmatrix
!  ! collect
!  ! extract
!
! SUBROUTINEs  :
!  fromQPtoCSR_real
!  fromQPtoCSR_cmplx
!  fromCSRtoQP_real
!  fromCSRtoQP_cmplx
!  
!  assembleCSRmatrix_real
!  assembleCSRmatrix_cmplx
!
!  countFromOne ( A , B )
!
!  sort_diff    (   ) 
!
!  collect_d    ( nVel , uu , pp , xx )
!  collect_z    ( nVel , uu , pp , xx )
!  extract_d    ( nVel , uu , pp , xx )
!  extract_z    ( nVel , uu , pp , xx )  
!  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

USE prep_mesh_p1p2_sp
USE global_variables


! INTERFACEs -----------------------------------------------------

INTERFACE fromQPtoCSR
  MODULE PROCEDURE fromQPtoCSR_real , &
                   fromQPtoCSR_cmplx
END INTERFACE fromQPtoCSR

INTERFACE fromCSRtoQP
  MODULE PROCEDURE fromCSRtoQP_real , &
                   fromCSRtoQP_cmplx
END INTERFACE fromCSRtoQP

INTERFACE assembleCSRmatrix
  MODULE PROCEDURE assembleCSRmatrix_real , &
                   assembleCSRmatrix_cmplx
END INTERFACE assembleCSRmatrix

! INTERFACE collect
!   MODULE PROCEDURE collect_d , &
!                    collect_z
! END INTERFACE collect
! 
! INTERFACE extract
!   MODULE PROCEDURE extract_d , &
!                    extract_z
! END INTERFACE extract
! 
! ----------------------------------------------------------------


CONTAINS

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE fromQPtoCSR_real(AAQP,AA)

  TYPE(CSR_MUMPS_MATRIX) , INTENT(IN) :: AAQP
  TYPE(CSR_MUMPS_MATRIX) :: AA
  
  WRITE(*,*) ' SIZE(AAQP%I) = ' , SIZE(AAQP%I)
  WRITE(*,*) ' SIZE(AAQP%J) = ' , SIZE(AAQP%J)
  ALLOCATE(AA%I      (SIZE(AAQP%I))      )
  ALLOCATE(AA%J      (SIZE(AAQP%J))      )
  ALLOCATE(AA%I_MUMPS(SIZE(AAQP%I_MUMPS)))
  ALLOCATE(AA%E      (SIZE(AAQP%E))      )
  
  AA%I       = AAQP%I - 1
  AA%I_MUMPS = AAQP%I_MUMPS - 1
  AA%J       = AAQP%J - 1
  AA%E       = AAQP%E
  
END SUBROUTINE fromQPtoCSR_real

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE fromQPtoCSR_cmplx(AAQP,AA)

  TYPE(CSR_MUMPS_COMPLEX_MATRIX) , INTENT(IN) :: AAQP
  TYPE(CSR_MUMPS_COMPLEX_MATRIX) :: AA
  
  ALLOCATE(AA%I      (SIZE(AAQP%I))      )
  ALLOCATE(AA%J      (SIZE(AAQP%J))      )
  ALLOCATE(AA%I_MUMPS(SIZE(AAQP%I_MUMPS)))
  ALLOCATE(AA%E      (SIZE(AAQP%E))      )
  
  AA%I       = AAQP%I - 1
  AA%I_MUMPS = AAQP%I_MUMPS - 1
  AA%J       = AAQP%J - 1
  AA%E       = AAQP%E
   
END SUBROUTINE fromQPtoCSR_cmplx

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE fromCSRtoQP_real(AA,AAQP)

  TYPE(CSR_MUMPS_MATRIX) , INTENT(IN) :: AA
  TYPE(CSR_MUMPS_MATRIX) :: AAQP
  
  ALLOCATE(AAQP%I      (SIZE(AA%I))      )
  ALLOCATE(AAQP%J      (SIZE(AA%I))      )
  ALLOCATE(AAQP%I_MUMPS(SIZE(AA%I_MUMPS)))
  ALLOCATE(AAQP%E      (SIZE(AA%E))      )
  
  AAQP%I       = AA%I + 1
  AAQP%I_MUMPS = AA%I_MUMPS + 1
  AAQP%J       = AA%J + 1
  AAQP%E       = AA%E
  
END SUBROUTINE fromCSRtoQP_real

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE fromCSRtoQP_cmplx(AA,AAQP)

  TYPE(CSR_MUMPS_COMPLEX_MATRIX) , INTENT(IN) :: AA
  TYPE(CSR_MUMPS_COMPLEX_MATRIX) :: AAQP
  
  ALLOCATE(AAQP%I      (SIZE(AA%I))      )
  ALLOCATE(AAQP%J      (SIZE(AA%J))      )
  ALLOCATE(AAQP%I_MUMPS(SIZE(AA%I_MUMPS)))
  ALLOCATE(AAQP%E      (SIZE(AA%E))      )
  
  AAQP%I       = AA%I + 1
  AAQP%I_MUMPS = AA%I_MUMPS + 1
  AAQP%J       = AA%J + 1
  AAQP%E       = AA%E
   
END SUBROUTINE fromCSRtoQP_cmplx

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE assembleCSRmatrix_real ( A , B , dRow , dCol , C )

! UP TO NOW, dRow and dCol MUST BE NON NEGATIVE
! INPUT AND OUTPUT MATRICES IN CSR FORMAT
  IMPLICIT NONE
  
  TYPE(CSR_MUMPS_MATRIX) , INTENT(IN)  :: A , B
  INTEGER                , INTENT(IN)  :: dRow , dCol
  
  TYPE(CSR_MUMPS_MATRIX) :: C
  INTEGER , DIMENSION(:) , ALLOCATABLE :: IBDEL , JBDEL
  INTEGER , DIMENSION(:) , ALLOCATABLE :: IA , JA , IB , JB
  INTEGER :: numColA
  
  INTEGER :: IND
  INTEGER :: I1
  
  IF ( (SIZE(A%E) .NE. SIZE(A%J)) .OR. (SIZE(A%E) .NE. SIZE(A%I_MUMPS)) ) THEN
     WRITE(*,*) 'INVALID FIRST CSR INPUT '
     STOP
  END IF
  IF ( (SIZE(A%E) .NE. SIZE(A%J)) .OR. (SIZE(A%E) .NE. SIZE(A%I_MUMPS)) ) THEN
     WRITE(*,*) 'INVALID FIRST CSR INPUT '
     STOP
  END IF
  
  ! DELTA ROW AND COL ( DELTA >= 0 )
  ALLOCATE( IBDEL (SIZE(B%I)+DROW) )
  IF ( dRow == 0 ) THEN
     IBDEL = B%I
  ELSE
     IBDEL(1:DROW) = 0
     IBDEL(DROW+1:SIZE(IBDEL)) = B%I
  END IF
  JBDEL = B%J + DCOL
  
  ! C%I
  ALLOCATE( C%I(MAX(SIZE(A%I),SIZE(IBDEL))) )
  ALLOCATE( IA(SIZE(C%I)) , IB(SIZE(C%I)) )
  IA = 0 ; IB = 0
  IA(1:SIZE(A%I  )) = A%I
  IB(1:SIZE(IBDEL)) = IBDEL
  IF ( SIZE( A%I ) .LT. SIZE(IA) ) THEN
     IA(SIZE( A%I )+1:SIZE(IA)) = IA   (SIZE(A%I)  )
  END IF
  IF ( SIZE(IBDEL) .LT. SIZE(IB) ) THEN
     IB(SIZE(IBDEL)+1:SIZE(IB)) = IBDEL(SIZE(IBDEL))
  END IF
  
  C%I = IA + IB

  ALLOCATE( C%J(SIZE(A%J)+SIZE(JBDEL)) , C%I_MUMPS(SIZE(A%J)+SIZE(JBDEL)) , C%E(SIZE(A%J)+SIZE(JBDEL)) )
  C%E = 0.0D0
  
  ! c%I_MUMPS
  DO I1 = 1 , SIZE(C%I) - 1
     C%I_MUMPS( 1+C%I(I1):C%I(I1+1) ) = I1 - 1
  END DO
  
  ! C%J and C%E
  numColA = MAXVAL(A%J) + 1
  DO I1 = 1 , SIZE(C%I) - 1
     
     C%J(1+C%I(I1):C%I(I1)+IA(I1+1)-IA(I1))   = A%J  (1+IA(I1):IA(I1+1)) 
     C%J(1+C%I(I1)+IA(I1+1)-IA(I1):C%I(I1+1)) = JBDEL(1+IB(I1):IB(I1+1)) + numColA
     
     C%E(1+C%I(I1):C%I(I1)+IA(I1+1)-IA(I1))   = A%E  (1+IA(I1):IA(I1+1)) 
     C%E(1+C%I(I1)+IA(I1+1)-IA(I1):C%I(I1+1)) = B%E  (1+IB(I1):IB(I1+1))
     
  END DO

END SUBROUTINE assembleCSRmatrix_real

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE assembleCSRmatrix_cmplx ( A , B , dRow , dCol , C )

! UP TO NOW, dRow and dCol MUST BE NON NEGATIVE
! INPUT AND OUTPUT MATRICES IN CSR FORMAT
  IMPLICIT NONE
  
  TYPE(CSR_MUMPS_COMPLEX_MATRIX) , INTENT(IN)  :: A , B
  INTEGER                , INTENT(IN)  :: dRow , dCol
  
  TYPE(CSR_MUMPS_COMPLEX_MATRIX) :: C
  INTEGER , DIMENSION(:) , ALLOCATABLE :: IBDEL , JBDEL
  INTEGER , DIMENSION(:) , ALLOCATABLE :: IA , JA , IB , JB
  INTEGER :: numColA
  
  INTEGER :: IND
  INTEGER :: I1
  
  IF ( (SIZE(A%E) .NE. SIZE(A%J)) .OR. (SIZE(A%E) .NE. SIZE(A%I_MUMPS)) ) THEN
     WRITE(*,*) 'INVALID FIRST CSR INPUT '
     STOP
  END IF
  IF ( (SIZE(A%E) .NE. SIZE(A%J)) .OR. (SIZE(A%E) .NE. SIZE(A%I_MUMPS)) ) THEN
     WRITE(*,*) 'INVALID FIRST CSR INPUT '
     STOP
  END IF
  
  ! DELTA ROW AND COL ( DELTA >= 0 )
  ALLOCATE( IBDEL (SIZE(B%I)+DROW) )
  IF ( dRow == 0 ) THEN
     IBDEL = B%I
  ELSE
     IBDEL(1:DROW) = 0
     IBDEL(DROW+1:SIZE(IBDEL)) = B%I
  END IF
  JBDEL = B%J + DCOL
  
  ! C%I
  ALLOCATE( C%I(MAX(SIZE(A%I),SIZE(IBDEL))) )
  ALLOCATE( IA(SIZE(C%I)) , IB(SIZE(C%I)) )
  IA = 0 ; IB = 0
  IA(1:SIZE(A%I  )) = A%I
  IB(1:SIZE(IBDEL)) = IBDEL
  IF ( SIZE( A%I ) .LT. SIZE(IA) ) THEN
     IA(SIZE( A%I )+1:SIZE(IA)) = IA   (SIZE(A%I)  )
  END IF
  IF ( SIZE(IBDEL) .LT. SIZE(IB) ) THEN
     IB(SIZE(IBDEL)+1:SIZE(IB)) = IBDEL(SIZE(IBDEL))
  END IF
  
  C%I = IA + IB

  ALLOCATE( C%J(SIZE(A%J)+SIZE(JBDEL)) , C%I_MUMPS(SIZE(A%J)+SIZE(JBDEL)) , C%E(SIZE(A%J)+SIZE(JBDEL)) )
  WRITE(*,*) ' SIZE(C%I_MUMPS) = ' , SIZE(C%I_MUMPS)
  C%E = 0.0D0
  
  ! c%I_MUMPS
  DO I1 = 1 , SIZE(C%I) - 1
     C%I_MUMPS( 1+C%I(I1):C%I(I1+1) ) = I1 - 1
  END DO
  
  ! C%J and C%E
  numColA = MAXVAL(A%J) + 1
  DO I1 = 1 , SIZE(C%I) - 1
     
    C%J(1+C%I(I1):C%I(I1)+IA(I1+1)-IA(I1))   = A%J  (1+IA(I1):IA(I1+1)) 
    C%J(1+C%I(I1)+IA(I1+1)-IA(I1):C%I(I1+1)) = JBDEL(1+IB(I1):IB(I1+1)) + numColA
    
    C%E(1+C%I(I1):C%I(I1)+IA(I1+1)-IA(I1))   = A%E  (1+IA(I1):IA(I1+1)) 
    C%E(1+C%I(I1)+IA(I1+1)-IA(I1):C%I(I1+1)) = B%E  (1+IB(I1):IB(I1+1))
     
  END DO

END SUBROUTINE assembleCSRmatrix_cmplx

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE countFromOne ( A , B )
     
     IMPLICIT NONE
     
     INTEGER , DIMENSION(:,:) , INTENT(IN)  :: A
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: B
     
     INTEGER , DIMENSION(:,:) , ALLOCATABLE :: ACOPY
     INTEGER :: I1 , I2 , I3
     INTEGER :: TOPV , MINV , IND , NUM
     
     IF ( ALLOCATED(B) ) THEN
        DEALLOCATE(B)
     END IF
     
     ALLOCATE(ACOPY(SIZE(A,1),SIZE(A,2)))
     ALLOCATE(B    (SIZE(A,1),SIZE(A,2)))
     B = 0
     ACOPY = A
     
     TOPV = MAXVAL(A) + 1
  !   MINV = MAXVAL(A) - 1
     IND  = 0
     NUM  = 0
     DO WHILE ( .NOT. (NUM == SIZE(ACOPY,1)*SIZE(ACOPY,2) ) )
        IND = IND + 1
        MINV = MINVAL(ACOPY)
        DO I2 = 1 , SIZE(ACOPY,1)
           DO I3 = 1 , SIZE(ACOPY,2)
              IF ( ACOPY(I2,I3) == MINV ) THEN
                 B(I2,I3)     = IND
                 ACOPY(I2,I3) = TOPV
                 NUM = NUM + 1
              ENDIF
           END DO
        END DO
     END DO
     
  
  END SUBROUTINE countFromOne 

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE sort_diff (a,  a_d, n_a_d)

!  sort in ascending order of the integer array  a  and generation
!  of the integer array  a_d  whose first  n_a_d  leading entries
!  contain different values in ascending order, while all the
!  remaining entries are set to zero

!  sorting by Shell's method.

   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(INOUT) :: a
   INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
   INTEGER,               INTENT(OUT)   :: n_a_d

   INTEGER :: n, na, inc, i, j, k, ia

   na = SIZE(a)

!  sort phase

   IF (na == 0) THEN
      n_a_d = 0
      RETURN
   ENDIF

   inc = 1
   DO WHILE (inc <= na)
      inc = inc * 3
      inc = inc + 1
   ENDDO

   DO WHILE (inc > 1)
      inc = inc/3
      DO i = inc + 1, na
         ia = a(i)
         j = i
         DO WHILE (a(j-inc) > ia)
            a(j) = a(j-inc)
            j = j - inc
            IF (j <= inc) EXIT
         ENDDO
         a(j) = ia
      ENDDO
   ENDDO

!  compression phase

   n = 1
   a_d(n) = a(1)
   DO k = 2, na
      IF (a(k) > a(k-1)) THEN
         n = n + 1
         a_d(n) = a(k)
      ENDIF
   ENDDO

   n_a_d = n

   a_d(n_a_d + 1 : na) = 0

END SUBROUTINE sort_diff

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE collect_d ( uu, pp,  xx )

! uu (velCmpnnts,np)   
! pp (np_L)            
! xx vector collecting the velocity components and the pressure on the
!    nodes  xx = (U1 , U2 , P)      , if  n_vel = 2
!           xx = (U1 , U2 , U3 , P) , if  n_vel = 3

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp
   
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: xx

   INTEGER :: n_vel
   INTEGER :: np, np_L
   INTEGER :: i1

   ! n. of components of the velocity (2 cart, 3 cyl axisym ...)
   n_vel = SIZE(uu,1)
    
   np   = SIZE(uu, 2)
   np_L = SIZE(pp)

   DO i1 = 1 , n_vel 
     xx(1+(i1-1)*np:i1*np) = uu(i1,:)
   ENDDO
   xx(n_vel*np+1:n_vel*np+np_L) = pp 

 
END SUBROUTINE collect_d

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE collect_z ( uu, pp,  xx )

! uu (velCmpnnts,np)   
! pp (np_L)            
! xx vector collecting the velocity components and the pressure on the
!    nodes  xx = (U1 , U2 , P)      , if  n_vel = 2
!           xx = (U1 , U2 , U3 , P) , if  n_vel = 3

   IMPLICIT NONE

   COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu
   COMPLEX(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp
   
   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: xx

   INTEGER :: n_vel
   INTEGER :: np, np_L
   INTEGER :: i1

   ! n. of components of the velocity (2 cart, 3 cyl axisym ...)
   n_vel = SIZE(uu,1)
    
   np   = SIZE(uu, 2)
   np_L = SIZE(pp)

   DO i1 = 1 , n_vel 
     xx(1+(i1-1)*np:i1*np) = uu(i1,:)
   ENDDO
   xx(n_vel*np+1:n_vel*np+np_L) = pp 

 
END SUBROUTINE collect_z

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE extract_d (xx,  uu,  pp)

   IMPLICIT NONE
    
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: xx
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: uu
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE, OPTIONAL :: pp
   
!  INTEGER :: np, np_L     PUBLIC in prep_mesh_p1p2_sp
   INTEGER :: n_vel
   INTEGER :: i1

   n_vel = velCmpnnts      ! <----- PUBLIC in global_variables

   IF ( ALLOCATED(uu) ) THEN ; DEALLOCATE(uu) ; ENDIF
   ALLOCATE(uu(n_vel,np)) 
  
   DO i1 = 1 , n_vel 
     uu(i1,:) = xx(1+(i1-1)*np:i1*np) 
   ENDDO
  
   IF (PRESENT(pp)) THEN 
      IF ( ALLOCATED(pp) ) THEN ; DEALLOCATE(pp) ; ENDIF
      ALLOCATE(pp(np_L))
      pp = xx(n_vel*np+1 : n_vel*np + np_L) 
   ENDIF
 
END SUBROUTINE extract_d

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE extract_z (xx,  uu,  pp)

   IMPLICIT NONE
    
   COMPLEX(KIND=8), DIMENSION(:),   INTENT(IN)  :: xx
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: uu
   COMPLEX(KIND=8), DIMENSION(:)  , ALLOCATABLE, OPTIONAL :: pp
   
   !INTEGER :: np, np_L
   INTEGER :: n_vel
   INTEGER :: i1

   n_vel = velCmpnnts      ! <----- PUBLIC in global_variables

   IF ( ALLOCATED(uu) ) THEN ; DEALLOCATE(uu) ; ENDIF
   ALLOCATE(uu(n_vel,np)) 
   DO i1 = 1 , n_vel 
     uu(i1,:) = xx(1+(i1-1)*np:i1*np) 
   ENDDO
  
   IF (PRESENT(pp)) THEN 
      IF ( ALLOCATED(pp) ) THEN ; DEALLOCATE(pp) ; ENDIF
      ALLOCATE(pp(np_L))
      pp = xx(n_vel*np+1 : n_vel*np + np_L) 
   ENDIF

 
END SUBROUTINE extract_z


END MODULE fem_miscellaneous
