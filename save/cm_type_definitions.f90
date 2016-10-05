MODULE cm_type_definitions

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module containing the TYPEs involved in the program
!
!
! last modified : 2016-03-15
!
! 2016-04-18 v1.0 : cylinder test ok
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! TYPEs:
!  Parameters
! TYPE parameters_type                                          UNUSED (!!)
!  Read Input and folder layout of the program
! TYPE  parfilenType 
! TYPE  folderType
!  Terms of the expansion for the Centre Manifold Reduction
! TYPE  orderType
! TYPE  powerExpansion
!  Influence terms from power expansion to nonlinear forcing terms
! TYPE  forcingCombinations_type
! TYPE  elementForcingInfluence_type
! TYPE  forcingInfluence_type
!
!
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IMPLICIT NONE

! -------------------------------------------------------------------------

  TYPE parameters_type
    
    INTEGER :: num
    REAL(KIND=8) , DIMENSION(:) , POINTER :: param
    
  END TYPE parameters_type

! -------------------------------------------------------------------------

  TYPE parFile
    
    CHARACTER(LEN=40) :: cm
    CHARACTER(LEN=40) :: lapack
    CHARACTER(LEN=40) :: mumps
    CHARACTER(LEN=40) :: arpack

  END TYPE parFile

! -------------------------------------------------------------------------

  TYPE baseSolFile
    
    CHARACTER(LEN=40) :: baseState
    CHARACTER(LEN=40) :: baseParam

  END TYPE baseSolFile

! -------------------------------------------------------------------------

  TYPE linearInputFile
    
    CHARACTER(LEN=40) :: massMatrix
    CHARACTER(LEN=40) :: jacobianMatrix
    CHARACTER(LEN=40) :: phi
    CHARACTER(LEN=40) :: psi
    CHARACTER(LEN=40) :: sigma

  END TYPE linearInputFile

! -------------------------------------------------------------------------

  TYPE nonlinearInputFile
    
    CHARACTER(LEN=40) :: abb

  END TYPE nonlinearInputFile

! -------------------------------------------------------------------------

  TYPE inputFile
    
    TYPE(parFile)            :: parfilen
    TYPE(baseSolFile)        :: base 
    TYPE(linearInputFile)    :: linear    
    TYPE(nonlinearInputFile) :: nonlinear    
    
  END TYPE inputFile

! -------------------------------------------------------------------------

  TYPE outputFile
    
    CHARACTER(LEN=40) :: Q
    CHARACTER(LEN=40) :: G
    CHARACTER(LEN=40) :: postpro

  END TYPE outputFile

! -------------------------------------------------------------------------
  
  TYPE meshFile
    
    CHARACTER(LEN=40) :: folder
    CHARACTER(LEN=40) :: mesh_name
    
  END TYPE meshFile

! -------------------------------------------------------------------------
  
  TYPE orderType
  
    INTEGER :: order
    INTEGER :: nterms
    INTEGER , DIMENSION(:,:) , POINTER :: indices 
  
  END TYPE orderType
  
! -------------------------------------------------------------------------
  
  TYPE powerExpansion
  
    INTEGER :: maxOrder
    INTEGER :: totTerms
    TYPE(orderType) , DIMENSION(:) , POINTER :: order
  
  END TYPE powerExpansion
  
! -------------------------------------------------------------------------

  TYPE forcingCombinations_type

   INTEGER , DIMENSION(:,:) , POINTER :: mat
   INTEGER :: Nt

  END TYPE forcingCombinations_type

! -------------------------------------------------------------------------

  TYPE elementForcingInfluence_type
    
    INTEGER , DIMENSION(:) , POINTER :: mind
    
  END TYPE elementForcingInfluence_type

! -------------------------------------------------------------------------

  TYPE forcingInfluence_type
    
    TYPE(elementForcingInfluence_type) , DIMENSION(:) , POINTER :: el
    
  END TYPE forcingInfluence_type

  
END MODULE cm_type_definitions
  

