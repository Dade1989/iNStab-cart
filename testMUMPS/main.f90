PROGRAM test

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! Solve Ax = b with MUMPS
! A can be singular; if the RHS is in Ran(A) the system
! has infinite solutions. Can MUMPS solve this problem?
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

  USE global_variables
  USE par_solve_mumps
  USE sparse_matrix_profiles
  
  IMPLICIT NONE

  INTEGER :: n , nnz
  TYPE(CSR_MUMPS_Matrix) :: A
  REAL(KIND=8) , DIMENSION(:) , ALLOCATABLE :: x , b

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

  ! A matrix : QP format = CSR starting from 1
!  !  A is singular
!  !  A =  1   0   1   0
!  !       2  -1   1   0 =  a1 | a2 | a3 | a4   with a3 = a1 + a2
!  !       0   1   1   0                             a4 = 0
!  !       0   0   0   0
! dim(Ker(A)) = 2
! Ker(A) = {  1    1   } = { k1 , k2 }
!          {  1    1   } 
!          { -1   -1   }
!          {  0    1   }
! A x = b 
!   1   0   1   0    x1     1     x = (1,2,0,0)+a1*k1+a2*k2
!   2  -1   1   0    x2  =  0
!   0   1   1   0    x3     2
!   0   0   0   0    x4     0

  n   = 5
  nnz = 8
  ALLOCATE(A%i      (n+1)) ; A%i       = (/ 1,3,6,8,9,9 /)
  ALLOCATE(A%i_mumps(nnz)) ; A%i_mumps = (/ 1,1,2,2,2,3,3,4 /)
  ALLOCATE(A%j      (nnz)) ; A%j       = (/ 1,3,1,2,3,2,3,4 /)
  ALLOCATE(A%e      (nnz)) ; A%e       = (/ 1.0,      1.0,    & 
                                            2.0,-1.0, 1.0,    &
                                                 1.0, 1.0,    &
                                                          0.0 /)
  
  ! b
  ALLOCATE(b(n)) ; b = (/ 1.0,0.0,2.0,0.0,0.0 /) 
  ALLOCATE(x(n)) ; x = b

!  ! A matrix : QP format = CSR starting from 1
!  !  A is singular
!  !  A =  1   0   1
!  !       2  -1   1  =  a1 | a2 | a3    with a3 = a1 + a2
!  !       0   1   1 
!  n   = 3;
!  nnz = 7;
!  ALLOCATE(A%i      (n+1)) ; A%i       = (/ 1,3,6,8 /)
!  ALLOCATE(A%i_mumps(nnz)) ; A%i_mumps = (/ 1,1,2,2,2,3,3 /)
!  ALLOCATE(A%j      (nnz)) ; A%j       = (/ 1,3,1,2,3,2,3 /)
!  ALLOCATE(A%e      (nnz)) ; A%e       = (/ 1.0,     1.0, &
!                                            2.0,-1.0,1.0, &
!                                                 1.0,1.0 /)
!  
!  ! b
!  ALLOCATE(b(n)) ; b = (/ 1.0,0.0,2.0 /) 
!  ALLOCATE(x(n)) ; x = b

! Calls to MUMPS
CALL par_mumps_master(INITIALIZATION ,1,A,0) 
CALL par_mumps_master(SYMBO_FACTOR   ,1,A,0) 
CALL par_mumps_master(NUMER_FACTOR   ,1,A,0)

! -- UPDATE MUMPS to version 5.0.0 (or 5.0.1) --
! for NULL PIVOT DETECTION in UNSYMMETRIC MATRICES
CALL par_mumps_master(DIRECT_SOLUTION,1,A,0,x) 
WRITE(*,*) " x = "
WRITE(*,*)   x

stop

   CALL par_mumps_master (FINALIZATION, 1, A, 0)

   ELSE
      DO WHILE ( 1 > 0 )
         CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
         CALL par_mumps (parMumpsJob, matrID)
      ENDDO
   ENDIF

CALL MPI_FINALIZE (MPI_COMM_WORLD, mpiIerr)

 
END PROGRAM test
