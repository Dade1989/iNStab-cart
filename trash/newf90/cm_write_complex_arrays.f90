MODULE cm_write_complex_arrays

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module containing the SUBROUTINEs to write varaible of type COMPLEX
!
! 2016-03-18
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! INTERFACE  :
!   writeComlex
! SUBROUTINEs:
!   writeComplex0 ( A , label , format )
!   writeComplex1 ( A , label , format , N )
!   writeComplex2 ( A , label , format , N )
!
! !!! change the input of the OUTPUT_FORMAT !!!
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  IMPLICIT NONE
  
  INTERFACE writeComplex
    
    MODULE PROCEDURE writeComplex0 , &
                     writeComplex1 , &
                     writeComplex2
    
  END INTERFACE writeComplex
  
  CONTAINS

! ----------------------------------------------------------------------

   SUBROUTINE writeComplex0( A , label , outForm )
! Inputs:
!  A       : complex vector to be printed on screen
!  outForm : format of the output
!  N       : number of columns on each printed row

  IMPLICIT NONE
  
  COMPLEX(KIND=8) , INTENT(IN) :: A
  CHARACTER(*) , INTENT(IN) :: outForm , label
  
  WRITE(*,'(A,A)',ADVANCE='no') label , " = " 
  IF ( AIMAG(A) .GE. 0.0d0 ) THEN
    WRITE(*,outForm)  & 
        DBLE(A) , " +" , ABS(AIMAG(A)) , "i  "
  ELSE
    WRITE(*,outForm)  & 
        DBLE(A) , " -" , ABS(AIMAG(A)) , "i  "
  ENDIF
  
  END SUBROUTINE writeComplex0
 
! ----------------------------------------------------------------------
  
  SUBROUTINE writeComplex1( A , label , outForm , N )
! Inputs:
!  A       : complex vector to be printed on screen
!  outForm : format of the output
!  N       : number of columns on each printed row

  IMPLICIT NONE
  
  COMPLEX(KIND=8) , DIMENSION(:) , INTENT(IN) :: A
  CHARACTER(*) , INTENT(IN) :: outForm , label
  INTEGER :: N , NEL
  
  INTEGER :: I1 , J , JEND
 
  WRITE(*,'(A,A)',ADVANCE='no') label , " = " 
  J = 0
  DO WHILE ( J .LT. SIZE(A,1) )
    JEND = MIN(J+N,SIZE(A,1))
    NEL = JEND-J
!    WRITE(*,'(A,I3,A,I3)',ADVANCE='NO') "Columns from " , J+1 , " to " , JEND
!    WRITE(*,*)
      DO I1 = 1 , NEL
        IF ( AIMAG(A(J+I1)) .GE. 0.0d0 ) THEN
          WRITE(*,outForm,ADVANCE='NO')  & 
              DBLE(A(J+I1)) , " +" , ABS(AIMAG(A(J+I1))) , "i  "
        ELSE
          WRITE(*,outForm,ADVANCE='NO')  & 
              DBLE(A(J+I1)) , " -" , ABS(AIMAG(A(J+I1))) , "i  "
        ENDIF
      ENDDO
!      WRITE(*,*) ""
    WRITE(*,*)
    J = J + N
  ENDDO
!  DO WHILE ( J .LT. SIZE(A,1) )
!    JEND = MIN(J+N,SIZE(A,1))
!    WRITE(*,*) " Columns from " , J+1 , " to " , JEND
!    WRITE(*,*) A(J+1:JEND)
!    WRITE(*,*)
!    J = J + N
!  ENDDO
  
  END SUBROUTINE writeComplex1
  
! ----------------------------------------------------------------------
  
  SUBROUTINE writeComplex2 ( A , label , outForm , N )
! Inputs:
!  A       : complex vector to be printed on screen
!  outForm : format of the output
!  N       : number of columns on each printed row

  IMPLICIT NONE
  
  COMPLEX(KIND=8) , DIMENSION(:,:) , INTENT(IN) :: A
  CHARACTER(LEN=*) , INTENT(IN) :: label , outForm
  INTEGER :: N , NEL
   
  INTEGER :: I1 , I2 ,  J , JEND

  WRITE(*,*)
  WRITE(*,*) label , " = "
  WRITE(*,*)

  J = 0
  DO WHILE ( J .LT. SIZE(A,2) )
    JEND = MIN(J+N,SIZE(A,2))
    NEL = JEND-J
    WRITE(*,'(A,I3,A,I3)',ADVANCE='NO') "Columns from " , J+1 , " to " , JEND
    WRITE(*,*)
    DO I1 = 1 , SIZE(A,1)
      DO I2 = 1 , NEL
        IF ( AIMAG(A(I1,J+I2)) .GE. 0.0d0 ) THEN
          WRITE(*,outForm,ADVANCE='NO')  & 
              DBLE(A(I1,J+I2)) , " +" , ABS(AIMAG(A(I1,J+I2))) , "i  "
        ELSE
          WRITE(*,outForm,ADVANCE='NO')  & 
              DBLE(A(I1,J+I2)) , " -" , ABS(AIMAG(A(I1,J+I2))) , "i  "
        ENDIF
      ENDDO
      WRITE(*,*) ""
    ENDDO
    WRITE(*,*)
    J = J + N
  ENDDO
  
  END SUBROUTINE writeComplex2
  
! ----------------------------------------------------------------------
  
END MODULE cm_write_complex_arrays
