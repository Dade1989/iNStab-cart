MODULE sparse_matrix_operations
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!

   USE global_variables ! MPI variables
   USE Dirichlet_Neumann
   USE prep_mesh_p1p2_sp

   IMPLICIT NONE


CONTAINS


!------------------------------------------------------------------------------
! subroutines for real matrices
!------------------------------------------------------------------------------
!
!
SUBROUTINE dAlinB_s (alpha, a, beta, b,  c)
!
!-----------------------------------------------------------------------
!         C  =  alpha A  +  beta B
!----------------------------------------------------------------------- 
! linearly combine two matrices
! Matrices are stored in compressed sparse row storage.
! The three matrices need to have the same sparsity pattern
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! a, b  = input real matrices in compressed sparse row format
!
! alpha,
! beta  = real coefficients of the linear combination
!
!
! on return:
!-----------
! c     = real matrix, containing the sum C  =  alpha A  +  beta B
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8),               INTENT(IN) :: alpha, beta
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a,     b
   ! output variables
   REAL(KIND=8), DIMENSION(:)             :: c
   ! local variables
   INTEGER :: i

!$OMP PARALLEL DO
   DO i = 1, SIZE(a)

      c(i) = alpha*a(i) + beta*b(i)

   END DO
!$OMP END PARALLEL DO

END SUBROUTINE dAlinB_s

!------------------------------------------------------------------------------

SUBROUTINE dAtimx (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array, containing the product y = A*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,      DIMENSION(:), INTENT(IN) :: ja, ia
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   REAL(KIND=8), DIMENSION(:)             :: y
   ! local variables
   REAL(KIND=8) :: yi
   INTEGER      :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
   IF ( SIZE(x) /= number_of_rows ) THEN
      WRITE(*,*) '--> dAtimx Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF
   ! A - y
   IF ( SIZE(y) /= number_of_rows ) THEN
      WRITE(*,*) '--> dAtimx Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF
   ! not necessary to check x-y consistency because A has to be square


   ! compute y <-- A*x
!$OMP PARALLEL DO
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x

      yi = 0d0

      DO p = ia(i), ia(i+1)-1

         yi = yi + a(p)*x(ja(p))

      END DO

      ! store result in y(i)

      y(i) = yi

   END DO
!$OMP END PARALLEL DO


END SUBROUTINE dAtimx

!------------------------------------------------------------------------------

SUBROUTINE dAtimx_T (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         transpose(A) times a vector
!----------------------------------------------------------------------- 
! multiplies the transpose of a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array, containing the product y = A'*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,      DIMENSION(:), INTENT(IN) :: ja, ia
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   REAL(KIND=8), DIMENSION(:)             :: y
   ! local variables
   INTEGER :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
   IF ( SIZE(x) /= number_of_rows ) THEN
      WRITE(*,*) '--> dAtimx_T Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF
   ! A - y
   IF ( SIZE(y) /= number_of_rows ) THEN
      WRITE(*,*) '--> dAtimx_T Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF
   ! not necessary to check x-y consistency because A has to be square

   ! compute y <-- A'*x
   y = 0

!  $OMP PARALLEL DO &    ! not working yet with big matrices,
!  $OMP PRIVATE(p)  &    ! don't know why
!  $OMP REDUCTION (+: y)
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x
      DO p = ia(i), ia(i+1)-1

         y(ja(p)) = y(ja(p)) + a(p)*x(i)

      END DO

   END DO
!  $OMP END PARALLEL DO


END SUBROUTINE dAtimx_T

!------------------------------------------------------------------------------

SUBROUTINE dEssM (a, ja, ia, n, b, jb, ib, i_mumpsb)
!
!-----------------------------------------------------------------------
!         Extract square sub-Matrix  nxn
!----------------------------------------------------------------------- 
! Extract a square submatrix B from matrix A
! Matrix A and B are stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! b, jb,
!    ib = output matrix in compressed sparse row format.
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,      DIMENSION(:), INTENT(IN) :: ja, ia
   INTEGER,                    INTENT(IN) :: n
   ! output variables
   REAL(KIND=8), DIMENSION(:), POINTER           :: b
   INTEGER,      DIMENSION(:), POINTER           :: jb, ib
   INTEGER,      DIMENSION(:), POINTER, OPTIONAL :: i_mumpsb
   ! local variables
   INTEGER      :: number_of_rows, i, p, c


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   !
   IF ( n > number_of_rows ) THEN
      WRITE(*,*) '--> dEssM Error: n > rows(A)' 
      WRITE(*,*) '    n = ', n, ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF


   ! count elements
   !
   c = 0
   DO i = 1, n

      DO p = ia(i), ia(i+1) - 1

         IF ( ja(p) <= n ) THEN

            c = c + 1

         ENDIF

      ENDDO

   ENDDO

   ALLOCATE( ib(n+1) )
   ALLOCATE( jb(c), b(c) )

   ! fill matrix
   !
   c = 0
   DO i = 1, n

      ib(i) = c + 1

      DO p = ia(i), ia(i+1) - 1

         IF ( ja(p) <= n ) THEN

            c = c + 1

            jb(c) = ja(p)
             b(c) =  a(p)

         ENDIF

      ENDDO

   ENDDO
   ib(n+1) = c + 1


   IF (PRESENT(i_mumpsb)) THEN

      ALLOCATE( i_mumpsb(SIZE(jb)) )

      DO i = 1, SIZE(ib) - 1
      
         DO p = ib(i), ib(i+1) - 1

            i_mumpsb(p) = i

         ENDDO

      END DO

   ENDIF

END SUBROUTINE dEssM

!------------------------------------------------------------------------------
!
SUBROUTINE dZeroM (A, n, B, i_mumpsb)
!
!-----------------------------------------------------------------------
TYPE(CSR_MUMPS_Matrix), INTENT(IN) :: A
INTEGER, INTENT(IN) :: n                  ! number of rows added at the bottom of A
TYPE(CSR_MUMPs_Matrix) :: B
INTEGER,      DIMENSION(:), POINTER, OPTIONAL :: i_mumpsb

ALLOCATE(B%e(SIZE(A%e)) , B%j(SIZE(A%j)))
B%e = A%e
B%j = A%j
ALLOCATE(B%i(SIZE(A%i)+n))
B%i(1:SIZE(A%i)) = A%i
B%i(SIZE(A%i)+1:SIZE(B%i)) = A%i(SIZE(A%i))


IF (PRESENT(i_mumpsb)) THEN
! why should I need i_mumpsb ??
END IF

END SUBROUTINE dZeroM

!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! subroutines for complex matrices
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ------------------------------------------------------------------------
SUBROUTINE orderSort ( a , ord )
   
   IMPLICIT NONE
   
   INTEGER , DIMENSION(:) , INTENT(IN) :: a
   INTEGER , DIMENSION(:) , ALLOCATABLE , INTENT(INOUT) :: ord
   
   INTEGER ,  DIMENSION(:) , ALLOCATABLE :: aCopy
   INTEGER :: M , MINOLD , EL
   INTEGER :: MINV , KMIN
   INTEGER :: I1 , i2
   
   ALLOCATE( ORD(SIZE(A)) , aCopy(SIZE(A)) )
   ORD = 0
   ACOPY = A
   
   M    = MAXVAL(A) + 1
   MINV = MINVAL(A) - 1
   EL = 0
   DO I1 = 1 , SIZE(A)
      
      MINOLD = MINV
      MINV = ACOPY(1)
      KMIN = 1
      DO I2 = 1 , SIZE(A)
         
         IF ( ACOPY(I2) .LE. MINV ) THEN
            MINV = ACOPY(I2)
            KMIN = I2
         END IF
         
      END DO
      IF ( MINV .NE. MINOLD ) THEN
         EL = EL + 1
      END IF
      ORD(KMIN)   = EL
      ACOPY(KMIN) = M
   END DO


END SUBROUTINE orderSort

! ------------------------------------------------------------------------------

SUBROUTINE sortDiff (a,  a_d, n_a_d)

!  sort in ascending order of the integer array  a  and generation
!  of the integer array  a_d  whose first  n_a_d  leading entries
!  contain different values in ascending order, while all the
!  remaining entries are set to zero

!  sorting by Shell's method.

   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(INOUT) :: a
   INTEGER, DIMENSION(:), ALLOCATABLE   :: a_d
   INTEGER :: n_a_d

   INTEGER :: n, na, inc, i, j, k, ia

   na = SIZE(a)
   ALLOCATE(a_d(na))
   
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

END SUBROUTINE sortDiff

! -----------------------------------------------------------------------------
  
   SUBROUTINE writeSparseMatrix ( AA )
   
   IMPLICIT NONE
   
   TYPE(CSR_MUMPS_COMPLEX_Matrix) :: AA
   INTEGER :: i1
   
   DO i1 = 1 , SIZE(AA%e)
   
      WRITE(*,*) AA%I_MUMPS(I1) , AA%J(I1) , AA%E(I1) 
   
   END DO
   
   END SUBROUTINE writeSparseMatrix
   
!------------------------------------------------------------------------------

   SUBROUTINE writeImumps_cmplx ( AA )
      
      IMPLICIT NONE
      
      TYPE(CSR_MUMPS_COMPLEX_Matrix) , INTENT(INOUT) :: AA
      INTEGER :: i1 , j1
      
      write(*,*) 'size(AA%i)       = ' , size(AA%i)
      write(*,*) 'size(AA%i_mumps) = ' , size(AA%i_mumps)
      
      DO i1 = 1, SIZE(AA%i)-1
         
         DO j1 = AA%i(i1), AA%i(i1+1)-1
            
            AA%i_mumps(j1) = i1 
!            write(*,*) j1+1
         ENDDO
         
      ENDDO
      
   END SUBROUTINE writeImumps_cmplx

! -------------------------------------------------------------------------------
!   SUBROUTINE hgInletVelocityB ( ndof , nblocks , nodesIn , nodesIn_Raw , A )
!   
!   IMPLICIT NONE 
!   
!   INTEGER, DIMENSION(:)  , INTENT(INOUT) :: nodesIn
!   INTEGER, DIMENSION(:)  , INTENT(INOUT) :: nodesIn_Raw
!   INTEGER, INTENT(IN) :: nblocks , ndof
!   TYPE(CSR_MUMPS_COMPLEX_Matrix) :: A
!   
!   INTEGER , DIMENSION(:) , ALLOCATABLE :: NODES
!   INTEGER , DIMENSION(:) , ALLOCATABLE :: ORD
!   INTEGER :: i1 , i2
!   integer :: n_a_d
!   
!   ! "EXTRA" INFLOW BOUNDARY MODIFICATIONS TO THE BBF IN ORDER TO CONSIDER THE RIGHT B.C.S ON THE AXIS (usually but not
!   !     always THE CENTRE OF THE INTAKE MANIFOLD
!   ALLOCATE( A%e(nblocks*SIZE(nodesIn_Raw)) ) ;  A%e = CMPLX(1.0d0,0.0d0,KIND=8)
!   ALLOCATE( A%i(nblocks*ndof + 1) )    ;
!   ALLOCATE( A%j(nblocks*SIZE(nodesIn_Raw)) ) ;
!   
!   CALL sortDiff  ( nodesIn_Raw , nodes , n_a_d)
!   CALL orderSort ( NODES , ord )!!!!!!!! THIS WAY     ORD = 1:size(nodes)
!   WRITE(*,*) ORD
!   
!   DO i1 = 1 , nblocks
!
!      DO i2 = 1 , SIZE(nodes)
!         
!         IF ( .NOT.(ANY( NODES(I2) == NODESIN ) ) ) THEN
!             A%e( (i1-1)*SIZE(nodes) + i2 ) = CMPLX(0.0d0,0.0d0,KIND=8)
!         END IF
!         
!         ! CORRECTION FOR THE CENTRE OF THE AXIS IF THE INFLOW MANIFOLD CONTAINS THE AXIS OF SYMMETRY
!         IF (BETA == 0) THEN
!           IF ( (NODES(i2) == INFLOW_CENTRE) .AND. (I1 .NE. 1) ) THEN
!             A%e( (i1-1)*SIZE(nodes) + i2 ) = CMPLX(0.0d0,0.0d0,KIND=8)
!           END IF
!         ELSEIF (BETA == 1) THEN
!           IF ( (NODES(i2) == INFLOW_CENTRE) .AND. (I1 == 1) ) THEN
!             A%e( (i1-1)*SIZE(nodes) + i2 ) = CMPLX(0.0d0,0.0d0,KIND=8)
!           END IF
!         ELSE
!           IF ( NODES(i2) == INFLOW_CENTRE ) THEN
!             A%e( (i1-1)*SIZE(nodes) + i2 ) = CMPLX(0.0d0,0.0d0,KIND=8)
!           END IF
!         END IF
!                                                                                            !!! matrix QP-style. Not CSR !!!
!         A%j( (i1-1)*SIZE(nodes) + i2 ) = (i1-1)*size(nodes) + ORD(i2)                      !!!!!!!!!!!!!!!!!!!! NO - 1  !!!
!         
!         A%i( (i1-1)*ndof+1 : (i1-1)*ndof+nodes(1) ) = (i1-1)*SIZE(nodes) + 1               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         IF ( i2 .LT. SIZE(nodes) ) THEN
!            A%i( (i1-1)*ndof + nodes(i2)+1 : (i1-1)*ndof + nodes(i2+1) ) = i2 + (i1-1)*SIZE(nodes) + 1
!         ELSE
!            A%i( (i1-1)*ndof + nodes(i2)+1 : i1*ndof ) = i2 + (i1-1)*SIZE(nodes) + 1
!         END IF
!         
!      END DO
!   
!   END DO
!   
!   A%i( nblocks*ndof + 1 ) = 3*SIZE(nodes) + 1
!   
!   END SUBROUTINE hgInletVelocityB 

! -------------------------------------------------------------------------------
   SUBROUTINE extractMatrix ( A , IROW , ICOL , B )      !!!! check the conventions used by QP and make this coherent !!!!!
   
   IMPLICIT NONE 
   
   TYPE(CSR_MUMPS_COMPLEX_Matrix) , INTENT(IN) :: A
   INTEGER, DIMENSION(:) , INTENT(IN) :: ICOL , IROW
   TYPE(CSR_MUMPS_COMPLEX_Matrix) , INTENT(INOUT) :: B
   
   INTEGER, DIMENSION(:) , ALLOCATABLE :: jBnew
!   INTEGER, DIMENSION(:) , ALLOCATABLE :: ICOL , IROW
   INTEGER :: naCol , naRow

   INTEGER :: p , q
   INTEGER :: nel
   COMPLEX(KIND=8) , DIMENSION(:) , ALLOCATABLE :: be
   INTEGER         , DIMENSION(:) , ALLOCATABLE :: bj , bjj , BjjOrd
!   COMPLEX(KIND=8), DIMENSION(:) , ALLOCATABLE :: 
!   INTEGER , DIMENSION(:), ALLOCATABLE :: AJord
   INTEGER :: lenRow
   
   INTEGER :: i1 , i2 
   

   nel = 0

!   CALL sortDiff ( IROWunsorted , IROW , naRow )
!   CALL sortDiff ( ICOLunsorted , ICOL , naCol )
!   write(*,*) ' naRow , naCol = ' , naCol , naRow
   
   p = SIZE(IROW) ; q = SIZE(ICOL)
   
   write(*,*) 'size(IROW) , SIZE(ICOL) = ' , SIZE(IROW) , SIZE(ICOL)
   
   ALLOCATE( B%i (p+1) )
   ALLOCATE( be(p*q) , bj(p*q) )
   be = cmplx(0.0d0,0.0d0,kind=8)
   bj = 0
!   WRITE(*,*) 'SIZE(A%J) = ' , SIZE(A%J)
!   CALL orderSort ( A%j , AJord )
   
   B%i(1) = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nel = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO i1 = 1 , p
      
      lenRow = A%i(IROW(i1)+1) - A%i(IROW(i1))                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! +1    e +0 ???
                                                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO i2 = 1 , lenRow                                                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
         IF ( ANY( ICOL == A%j( A%i(IROW(i1))+i2 -1) ) ) THEN
            
            nel = nel + 1
            be(nel) = A%e( A%i(IROW(i1))+i2 -1 )
            bj(nel) = A%j( A%i(IROW(i1))+i2 -1 ) - MINVAL(A%j) +1  !(1)
            
         END IF
         
      END DO
      
      B%i(i1+1) = nel + 1
      
   END DO
   
   ALLOCATE( B%e(nel) ) ; B%e = be(1:nel)
   ALLOCATE( Bjj(nel) ) ; Bjj = bj(1:nel)
   ALLOCATE( B%i_mumps(nel) ) ; B%i_mumps = 1
   
   WRITE(*,*) ' NEL = ' ,  NEL
   
   WRITE(*,*) 'Starting sorting ...'
   CALL orderSort ( Bjj , BjjOrd )
      
   WRITE(*,*) '   sorting done !'
   write(*,*) 'size(BJJORD)' , SIZE(BJJORD)
   ALLOCATE( B%j(nel) ) ; B%j = bjjord             !!!!!!!!!!!!!!!!

   CALL  writeImumps_cmplx ( B )
  
   END SUBROUTINE extractMatrix  
      
! -------------------------------------------------------------------------------

!!SUBROUTINE sparseMatMul( A , B , C )                                                   !!! TOO SLOW !!!
!!   
!!   TYPE(CSR_MUMPS_COMPLEX_Matrix) , INTENT(IN) :: A , B
!!   TYPE(CSR_MUMPS_COMPLEX_Matrix)              :: C
!!   
!!   COMPLEX(KIND=8) , DIMENSION(:) , ALLOCATABLE :: Ce , Cj
!!   
!!   TYPE(CSR_MUMPS_COMPLEX_Matrix) :: BT
!!   TYPE(CSR_MUMPS_COMPLEX_Vector) :: av , bv
!!   
!!   COMPLEX(KIND=8) :: z
!!   
!!   INTEGER :: ABEG , AEND , BTBEG , BTEND
!!   INTEGER :: i1 , i2
!!   INTEGER :: numElC
!!   INTEGER :: ind
!!   
!!   CALL sparseTranspose ( B , BT )
!!   
!!   ALLOCATE( Ce ( (SIZE(A%I)-1) * (SIZE(BT%I)-1) ) ) ; CE  = 0.0D0 
!!   ALLOCATE( C%I(  SIZE(A%I) ) )                     ; C%I = 0
!!   ALLOCATE( Cj ( (SIZE(A%I)-1) * (SIZE(BT%I)-1) ) ) ; CJ  = 0
!!   
!!   numElC = 0
!!   C%i(1) = 0
!!   
!!   DO i1 = 1 , SIZE(A%I)-1            ! ROWS OF A
!!      C%i(i1+1) = C%i(i1)
!!!      IND = 0
!!      ABEG = A%J( A%I(I1)+1 ) +1
!!      AEND = A%J( A%I(I1+1) ) +1
!!!      write(*,*) i1 , ABEG , AEND      
!!      IF ( AEND .GE. ABEG ) THEN
!!         ALLOCATE( av%e (AEND-ABEG+1) )   ;  av%e = A%E( A%I(I1)+1 : A%I(I1+1) )
!!         ALLOCATE( av%j (AEND-ABEG+1) )   ;  av%j = A%J( A%I(I1)+1 : A%I(I1+1) )
!!!         write(*,*) ' i1 = ' , i1 , ' av%e ='
!!!         write(*,*) av%e
!!         DO i2 = 1 , SIZE(BT%I)-1        ! ROWS OF BT = COLS OF B
!!            
!!            BTBEG = BT%J(BT%I(I2)+1)+1
!!            BTEND = BT%J(BT%I(I2+1))+1
!!            IF ( BTEND .GE. BTBEG ) THEN
!!               IND = 0
!!               ALLOCATE( bv%e (BTEND-BTBEG+1) )   ;  bv%e = BT%E( BT%I(I2)+1 : BT%I(I2+1) )
!!               ALLOCATE( bv%j (BTEND-BTBEG+1) )   ;  bv%j = BT%J( BT%I(I2)+1 : BT%I(I2+1) )
!!               
!!               CALL sparseScalarProduct ( av , bv , z , IND )
!!!               write(*,*) ind
!!!               WRITE(*,*) ' i1 , i2 , z = ' , i1 , i2 , z
!!               
!!               IF ( IND .NE. 0 ) THEN
!!                  numElC = numElC + 1
!!                  Ce(numElC) = z
!!                  C%i(i1+1)  = C%i(i1+1) + 1
!!                  Cj(numElC) = i2 - 1
!!               END IF
!!               
!!               DEALLOCATE(bv%e,bv%j)
!!            END IF      
!!         END DO
!!         DEALLOCATE(av%e,av%j) 
!!      ELSE
!!         
!!      END IF
!!      
!!      
!!      
!!   END DO
!!   
!!   ALLOCATE( C%E ( numElC ) )   ; C%E = Ce(1:numElC)
!!   ALLOCATE( C%J ( numElC ) )   ; C%J = Cj(1:numElC)
!!   
!!   
!!   
!!END SUBROUTINE sparseMatMul

!!!------------------------------------------------------------------------------

!!SUBROUTINE sparseScalarProduct( a , b , z , ind )                                 !!! TOO SLOW !!!
!!!
!!! IN  : A
!!! OUT : B = A^T

!!   IMPLICIT NONE
!!   
!!   TYPE(CSR_MUMPS_COMPLEX_Vector) , INTENT(IN) :: a , b
!!   COMPLEX(KIND = 8) :: z
!!   INTEGER :: ind
!!   
!!   INTEGER :: i1 , i2
!!   
!!   z = 0.0d0
!!   DO i1 = 1 , SIZE(a%j)
!!      DO i2 = 1 , SIZE(b%j)
!!         IF (a%j(i1) == b%j(i2)) THEN
!!            z = z + a%e(i1) * b%e(i2)
!!         ENDIF
!!      END DO
!!   END DO
!!   
!!   IF ( z .NE. 0.0d0 ) THEN
!!      ind = 1
!!   END IF
!!   
!!   END SUBROUTINE sparseScalarProduct


!------------------------------------------------------------------------------
!
! -
SUBROUTINE sparseTranspose ( A , B )       !!!!!!!!!!!!!!!!!!!! check the conventions used by QP and make this coherent !!!!!!!
!
! IN  : A
! OUT : B = A^T

  IMPLICIT NONE
   
   TYPE(CSR_MUMPS_COMPLEX_Matrix) , INTENT(IN) :: A
   TYPE(CSR_MUMPS_COMPLEX_Matrix)              :: B
   INTEGER , DIMENSION(:) , ALLOCATABLE :: vA                       ! VCTR CONTAINING THE COLUMN INDICES OF A
   
   INTEGER , DIMENSION(:) , ALLOCATABLE :: setk
   INTEGER :: numk
   INTEGER :: nvA , nel
   INTEGER :: ntot
   INTEGER :: col , rowN
   
   INTEGER :: i1 , i2 , i3
   
   ntot = 0
   
   ! Build va
   ALLOCATE( vA(SIZE(A%e)) ) 
   vA = 0
   nvA = 0
   DO i1 = 1 , SIZE(A%i) -1
      nel = A%i(i1+1) - A%i(i1)
      vA(nvA+1:nvA+nel) = i1 -1
      nvA = nvA + nel
   END DO
   
!   WRITE(*,*) ' vA = '
!   WRITE(*,*)   vA
   
   !
   ALLOCATE( B%e (SIZE(A%e)) )      ; B%e = 0.0d0
   ALLOCATE( B%i (MAXVAL(A%j)+2) )  ; B%i = 0
   ALLOCATE( B%j (SIZE(A%e)) )      ; B%j = 0
   
   ALLOCATE( setk(SIZE(A%j)) )
   
   B%i(1) = 0
   DO i1 = 0 , MAXVAL(A%j)
      
      rowN = i1 + 1
      
      ! Find k s.t. A%j(k) = i1
      setk = 0
      numk = 0
      DO i2 = 1 , SIZE(A%j)
         IF ( A%j(i2) == i1 ) THEN
            numk = numk + 1
            setk(numk) = i2
         ENDIF
      END DO
      ! setk = ...
      ! numk = ...
      B%e( ntot+1 : ntot+numk ) = A%e( setk(1:numk) )
      B%i( i1+2 ) = ntot + numk
      B%j( ntot+1 : ntot+numk ) = vA ( setk(1:numk) )
      ntot = ntot + numk
      
   END DO
   
   DEALLOCATE( setk )
   
   END SUBROUTINE sparseTranspose

! ----------------------------------------------------------------------

SUBROUTINE zAlinB_s (alpha, a, beta, b,  c)
!
!-----------------------------------------------------------------------
!         C  =  alpha A  +  beta B
!----------------------------------------------------------------------- 
! linearly combine two matrices
! Matrices are stored in compressed sparse row storage.
! The three matrices need to have the same sparsity pattern
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! a, b  = input complex matrices in compressed sparse row format
!
! alpha,
! beta  = complex coefficients of the linear combination
!
!
! on return:
!-----------
! c     = complex matrix, containing the sum C  =  alpha A  +  beta B
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8),               INTENT(IN) :: alpha, beta
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a,     b
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:)             :: c
   ! local variables
   INTEGER :: i

!$OMP PARALLEL DO
   DO i = 1, SIZE(a)

      c(i) = alpha*a(i) + beta*b(i)

   END DO
!$OMP END PARALLEL DO

END SUBROUTINE zAlinB_s

!------------------------------------------------------------------------------

SUBROUTINE zAtimx (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! on entry:
!----------
! x     = complex array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = complex array, containing the product y = A*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,         DIMENSION(:), INTENT(IN) :: ja, ia
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:)             :: y
   ! local variables
   COMPLEX(KIND=8) :: yi
   INTEGER         :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
!!!!!!!!!!!!!!!!   IF ( SIZE(x) /= number_of_rows ) THEN
!!!          !!!      WRITE(*,*) '--> zAtimx Error: wrong A - x size' 
!!!  FALSO   !!!      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
!!!          !!!      WRITE(*,*) '    STOP.'
!!!!!!!!!!!!!!!!      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
!!!!!!!!!!!!!!!!   END IF
   ! A - x (if A is CSR format, the only meaningful check is )
   IF ( SIZE(x) < MAXVAL(ja)) THEN ! ERROR
      WRITE(*,*) '--> zAtimx Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' min(SIZE(A,2)) = ', MAXVAL(ja)
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ELSEIF ( SIZE(x) > MAXVAL(ja)) THEN
!      WRITE(*,*) '--> zAtimx Warning: wrong A - x size' 
!      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' min(SIZE(A,2)) = ', MAXVAL(ja)
!      WRITE(*,*) '    The operation is allowed, but some entries in x are ignored.'
   END IF
   ! A - y
   IF ( SIZE(y) /= number_of_rows ) THEN
      WRITE(*,*) '--> zAtimx Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF
   ! not necessary to check x-y consistency because A has to be square

   ! compute y <-- A*x
!$OMP PARALLEL DO
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x

      yi = CMPLX(0d0, 0d0, KIND=8)

      DO p = ia(i), ia(i+1)-1

         yi = yi + a(p)*x(ja(p))

      END DO

      ! store result in y(i)

      y(i) = yi

   END DO
!$OMP END PARALLEL DO


END SUBROUTINE zAtimx

!------------------------------------------------------------------------------

SUBROUTINE zAtimx_T (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         transpose(A) times a vector
!----------------------------------------------------------------------- 
! on entry:
!----------
! x     = complex array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = complex array, containing the product y = A'*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,         DIMENSION(:), INTENT(IN) :: ja, ia
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:)             :: y
   ! local variables
   INTEGER :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
   IF ( SIZE(x) /= number_of_rows ) THEN
      WRITE(*,*) '--> zAtimx_T Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF
   ! A - y
      ! check dimension consistency
!!!!!!!!!!!!!   IF ( SIZE(y) /= number_of_rows ) THEN
!!!       !!!      WRITE(*,*) '--> zAtimx_T Error: wrong A - y size' 
!!! FALSO !!!      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
!!!       !!!      WRITE(*,*) '    STOP.'
!!!!!!!!!!!!!      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
!!!!!!!!!!!!!   END IF
   ! A - y (if A is CSR format, the only meaningful check is )
   IF ( SIZE(y) < MAXVAL(ja)) THEN ! ERROR
      WRITE(*,*) '--> zAtimx Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' min(SIZE(A,2)) = ', MAXVAL(ja)
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ELSEIF ( SIZE(y) > MAXVAL(ja)) THEN
!      WRITE(*,*) '--> zAtimx Warning: wrong A - y size' 
!      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' min(SIZE(A,2)) = ', MAXVAL(ja)
!      WRITE(*,*) '    The operation is allowed, but some entries in y are ignored.'
   END IF

   ! compute y <-- A'*x
   y = CMPLX(0d0, 0d0, KIND=8)
   
!  $OMP PARALLEL DO &    ! not working yet with big matrices,
!  $OMP PRIVATE(p)  &    ! don't know why
!  $OMP REDUCTION (+: y)
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x
      DO p = ia(i), ia(i+1)-1

         y(ja(p)) = y(ja(p)) + a(p)*x(i)

      END DO

   END DO
!  $OMP END PARALLEL DO


END SUBROUTINE zAtimx_T

!------------------------------------------------------------------------------

SUBROUTINE zAtimx_H (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         hermitian(A) times a vector
!----------------------------------------------------------------------- 
! on entry:
!----------
! x     = complex array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = complex array, containing the product y = A^H*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,         DIMENSION(:), INTENT(IN) :: ja, ia
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:)             :: y
   ! local variables
   INTEGER :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
   IF ( SIZE(x) /= number_of_rows ) THEN
      WRITE(*,*) '--> zAtimx_H Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF
   ! A - y
      ! check dimension consistency
!!!!!!!!!!!!!   IF ( SIZE(y) /= number_of_rows ) THEN
!!!       !!!      WRITE(*,*) '--> zAtimx_H Error: wrong A - y size' 
!!! FALSO !!!      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
!!!       !!!      WRITE(*,*) '    STOP.'
!!!!!!!!!!!!!      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
!!!!!!!!!!!!!   END IF
   ! A - y (if A is CSR format, the only meaningful check is )
   IF ( SIZE(y) < MAXVAL(ja)) THEN ! ERROR
      WRITE(*,*) '--> zAtimx Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' min(SIZE(A,2)) = ', MAXVAL(ja)
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ELSEIF ( SIZE(y) > MAXVAL(ja)) THEN
!      WRITE(*,*) '--> zAtimx Warning: wrong A - y size' 
!      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' min(SIZE(A,2)) = ', MAXVAL(ja)
!      WRITE(*,*) '    The operation is allowed, but some entries in y are ignored.'
   END IF

   ! compute y <-- A'*x
   y = CMPLX(0d0, 0d0, KIND=8)
   
!  $OMP PARALLEL DO &    ! not working yet with big matrices,
!  $OMP PRIVATE(p)  &    ! don't know why
!  $OMP REDUCTION (+: y)
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x
      DO p = ia(i), ia(i+1)-1

         y(ja(p)) = y(ja(p)) + DCONJG(a(p))*x(i)

      END DO

   END DO
!  $OMP END PARALLEL DO


END SUBROUTINE zAtimx_H

!------------------------------------------------------------------------------

SUBROUTINE zEssM (a, ja, ia, n, b, jb, ib, i_mumpsb)
!
!-----------------------------------------------------------------------
!         Extract square sub-Matrix  nxn
!----------------------------------------------------------------------- 
! Extract a square submatrix B from matrix A
! Matrix A and B are stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 16/03/2014
!
! on entry:
!----------
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! b, jb,
!    ib = output matrix in compressed sparse row format.
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,         DIMENSION(:), INTENT(IN) :: ja, ia
   INTEGER,                       INTENT(IN) :: n
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:), POINTER        :: b
   INTEGER,      DIMENSION(:), POINTER           :: jb, ib
   INTEGER,      DIMENSION(:), POINTER, OPTIONAL :: i_mumpsb
   ! local variables
   INTEGER      :: number_of_rows, i, p, c


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   !
   IF ( n > number_of_rows ) THEN
      WRITE(*,*) '--> zEssM Error: n > rows(A)' 
      WRITE(*,*) '    n = ', n, ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF


   ! count elements
   !
   c = 0
   DO i = 1, n

      DO p = ia(i), ia(i+1) - 1

         IF ( ja(p) <= n ) THEN

            c = c + 1

         ENDIF

      ENDDO

   ENDDO

   ALLOCATE( ib(n+1) )
   ALLOCATE( jb(c), b(c) )

   ! fill matrix
   !
   c = 0
   DO i = 1, n

      ib(i) = c + 1

      DO p = ia(i), ia(i+1) - 1

         IF ( ja(p) <= n ) THEN

            c = c + 1

            jb(c) = ja(p)
             b(c) =  a(p)

         ENDIF

      ENDDO

   ENDDO
   ib(n+1) = c + 1


   IF (PRESENT(i_mumpsb)) THEN

      ALLOCATE( i_mumpsb(SIZE(jb)) )

      DO i = 1, SIZE(ib) - 1
      
         DO p = ib(i), ib(i+1) - 1

            i_mumpsb(p) = i

         ENDDO

      END DO

   ENDIF

END SUBROUTINE zEssM
!
!

!------------------------------------------------------------------------------
!
SUBROUTINE zZeroM (A, n, B, i_mumpsb)
!
!-----------------------------------------------------------------------
TYPE(CSR_MUMPS_COMPLEX_Matrix), INTENT(IN) :: A
INTEGER, INTENT(IN) :: n                  ! number of rows added at the bottom of A
TYPE(CSR_MUMPS_COMPLEX_Matrix) :: B
INTEGER,      DIMENSION(:), POINTER, OPTIONAL :: i_mumpsb

ALLOCATE(B%e(SIZE(A%e)) , B%j(SIZE(A%j)))
B%e = A%e
B%j = A%j
ALLOCATE(B%i(SIZE(A%i)+n))
B%i(1:SIZE(A%i)) = A%i
B%i(SIZE(A%i)+1:SIZE(B%i)) = A%i(SIZE(A%i))


IF (PRESENT(i_mumpsb)) THEN
! why should I need i_mumpsb ??
END IF

END SUBROUTINE zZeroM

!==============================================================================

END MODULE sparse_matrix_operations
