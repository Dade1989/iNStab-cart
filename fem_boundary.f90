MODULE fem_boundary

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Module containing the ROUTINEs to build boundary and sides structures
!
! 2016-05-18
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! TYPEs :
!  sidesId_str                                  in global_variables.f90
!
! SUBROUTINEs :
!  ... read_and_apply_boundary_conditions            TO NOW in main.f90
!  fromJJStoJJSw          ( jjs , mmsw , jjswt , jjsw )
!  boundaryElementsOfSides( iMms , iSides , dSides , dMms )
!  unionFromMatToVct      ( ) 
!  inflowNodesOnly        ( )
!  innerNodesOnly         ( )
!  build_Id_sides_str     ( id_sides , sides_l_str )
!  build_logic_Wal_In     ( velCmpnnts , n.of sides , sides_l_str,
!                                                    logicWa, logicIn )
!  boundaryElementsOfSides(iMms,iSides,dSides,dMms)
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

USE global_variables
USE fem_miscellaneous

IMPLICIT NONE


! PUBLIC definitions
! ...

CONTAINS

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE fromJJStoJJSw ( JJS , MMSW , JJSWT , JJSW )
! build the boundary numbering

   IMPLICIT NONE
   
   INTEGER , DIMENSION(:,:) , INTENT(IN)  :: JJS
   INTEGER , DIMENSION(:)   , INTENT(IN)  :: MMSW
   INTEGER , DIMENSION(:,:) , ALLOCATABLE :: JJSWT , JJSW
   
   IF ( ALLOCATED(JJSWT) ) THEN
      DEALLOCATE(JJSWT)
   END IF
   IF ( ALLOCATED(JJSW) ) THEN
      DEALLOCATE(JJSW)
   END IF
   
   ALLOCATE( JJSWT(SIZE(JJS,1),SIZE(MMSW)) , JJSW(SIZE(JJS,1),SIZE(MMSW)) )
   JJSWT = 0 ; JJSW = 0
   JJSWT = JJS(:,MMSW)
   
   CALL countFromOne ( JJSWT , JJSW )
   
END SUBROUTINE 

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE  boundaryElementsOfSides(iMms,iSides,dSides,dMms)

! Routine that selects the boundary elements belonging to the desired
!  sides from the complete set
! INPUTs  :
!  iMms   = array of the numbering of all the boundary elements 
!         =   1:n.boundary elements               (n.boundary elements)
!  iSides = array of the side label of all the boundary elements :  
!             entries from 1 to n.sides           (n.boundary elements)
!  dSides = array of the desired labels : possible entries from 1 to n.
!           sides                                     (n.desired sides)
! OUTPUTs :
!  dMms   = array of the numbering of the desired boundary elements
!         = entries from 1 to n.boundary elements (n.des.boundary elem)
!
   IMPLICIT NONE
   
   INTEGER, DIMENSION(:), INTENT(IN)  :: iMms , iSides , dSides
   INTEGER, DIMENSION(:), ALLOCATABLE :: dMms
   INTEGER :: nDesElements , kDes
   
   INTEGER :: i1 , i2
   
   ! Find the dimension of the array dMms
   nDesElements = 0
   
   DO i1 = 1 , SIZE(iMms)
   
      DO i2 = 1 , SIZE(dSides)
        
         IF ( iSides(i1) == dSides(i2) ) THEN
            
            nDesElements = nDesElements + 1
            
         END IF
         
      END DO
   
   END DO
   
   ! Allocate dMms and fill it
   IF ( ALLOCATED(dMms) ) THEN
      DEALLOCATE(dMms)
   END IF
   
   ALLOCATE(dMms(nDesElements))
   kDes = 0
   DO i1 = 1 , SIZE(iMms)
   
      DO i2 = 1 , SIZE(dSides)
        
         IF ( iSides(i1) == dSides(i2) ) THEN
            
            kDes = kDes + 1
            
            dMms(kDes) = iMms(i1)
            
         END IF
         
      END DO
   
   END DO
   
   
   
END SUBROUTINE boundaryElementsOfSides

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE unionFromMatToVct(matrix,union)
! Routine that makes the set union of the elements collected in a matrix
! INPUT and OUTPUT are INTEGER greater than 0
! IN  : matrix
! OUT : union


   IMPLICIT NONE

   INTEGER, DIMENSION(:,:), INTENT(IN)  :: matrix
   INTEGER, DIMENSION(:)  , ALLOCATABLE , INTENT(INOUT) :: union
   
   INTEGER, DIMENSION(:)  , ALLOCATABLE :: unionUns
   INTEGER :: NAD
   INTEGER :: m , n
   INTEGER :: mn , nUnion
   INTEGER, DIMENSION(:)  , ALLOCATABLE :: unionExtended
   INTEGER :: ioAlreadyExist
   
   INTEGER :: i1 , i2 , i3
   
   
   m = SIZE(matrix,1) ;   n = SIZE(matrix,2) 
   mn = m * n
   
   ALLOCATE(unionExtended(mn))
   unionExtended = 0
   nUnion = 0
   DO i1 = 1 , m
      DO i2 = 1 , n
         ioAlreadyExist = 0
         DO i3 = 1 , nUnion                  ! " Incremental Search"
            IF ( unionExtended(i3) ==  matrix(i1,i2)) THEN
               ioAlreadyExist = 1
            END IF
         END DO
         IF ( ioAlreadyExist == 0 ) THEN
            nUnion = nUnion + 1
            unionExtended(nUnion) = matrix(i1,i2)
         ENDIF
      END DO
   END DO
   
   ALLOCATE(union(nUnion),unionUns(nUnion))
   unionUns = unionExtended(1:nUnion)
   
   CALL sort_diff(unionUns,union,NAD)
   
   
   DEALLOCATE(unionExtended)
   
   
   
END SUBROUTINE unionFromMatToVct

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE inflowNodesOnly(JJS_ALL,MMS_SID,NOD_RAW,NOD_AXIS,NODES, &
                                   IN_CENTRE)   !,nodeInflowCentre)

   IMPLICIT NONE
! IN  :   
!    JJS_ALL :: ARRAY CONTAINING THE GLOBAL NUMBERING OF ALL THE
!               BOUNDARY NODES ; JJS(NDOF_1EL , NSIDES)
!    MMS_SID :: ID OF THE BOUNDARY ELEMENTS DESIRED
!       EX : JJS_ALL(:,MMS_SID) = GLOBAL NUMBERING OF THE NODES ON THE
!                                  DESIRED SIDES
!    NOD_RAW ::  LIST TO BE CORRECTED
!    
! OUT :
!    NODES :: CORRECTED LIST

   INTEGER, DIMENSION(:,:) , INTENT(IN) :: JJS_ALL
   INTEGER, DIMENSION(:)   , INTENT(IN) :: MMS_SID
   INTEGER, DIMENSION(:)   , INTENT(IN) :: NOD_AXIS
   INTEGER, DIMENSION(:)   , INTENT(IN) :: NOD_RAW
   
   INTEGER, DIMENSION(:)   , ALLOCATABLE :: NODES
   INTEGER :: IN_CENTRE
   
   INTEGER, DIMENSION(:,:) , ALLOCATABLE :: JJS_SWAP
   INTEGER, DIMENSION(:)   , ALLOCATABLE :: NODES_L
   INTEGER :: N   ! NUMBER OF NODES TO BE REMOVED
   INTEGER :: I1

   ! INFLOW_CENTRE : ID OF THE NODES OF THE INTAKE MANIFOLD BELONGING
   !   TO THE AXIS OF SYMMETRY
   ! INITIALIZED TO 0: IF THE AXIS DOESN'T CROSS THE INTAKE MANIFOLD,
   !   INFLOW_CENTRE KEEPS 0 VALUE
   IN_CENTRE = 0
   WRITE(*,*) 'IN_CENTRE = ' , IN_CENTRE
   
   ALLOCATE(NODES_L(SIZE(NOD_RAW))) ; NODES_L = 0
   ALLOCATE(JJS_SWAP(SIZE(JJS_ALL,1) , SIZE(JJS_ALL,2)))
   JJS_SWAP = JJS_ALL
   JJS_SWAP(:,MMS_SID) = 0
   N = 0
   DO I1 = 1 , SIZE(NOD_RAW)
      
      IF ( ( .NOT.( ANY( NOD_RAW(I1) == JJS_SWAP(:,:) ) ) ) .OR. &
           (      ( ANY( NOD_RAW(I1) == NOD_AXIS(:)   ) ) )  ) THEN
         IF ( ANY( NOD_RAW(I1) == NOD_AXIS(:) ) ) THEN
            IN_CENTRE = NOD_RAW(I1)
            WRITE(*,*) "INFLOW_CENTRE = " , IN_CENTRE
         END IF
         N = N + 1
         NODES_L(N) = NOD_RAW(I1)
      END IF
      
   END DO 
   WRITE(*,*) 'IN_CENTRE = ' , IN_CENTRE
   
   ALLOCATE(NODES(N))
   NODES = NODES_L(1:N)

END SUBROUTINE inflowNodesOnly

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE   build_Id_sides_str ( id_Sides , sides_l_str )

! INPUT  :
! . id_Sides (vector of n.sides char):
!     sides label : 'i' inflow , 'o' outflow , 'w' wall , 'a' axis
! OUTPUT :
! . sides_l_str
!     

  IMPLICIT NONE
  
  CHARACTER(LEN=1) , DIMENSION(:) , INTENT(IN) :: id_Sides
  TYPE(sidesId_str) :: sides_l_str
  
  INTEGER :: numIn , numWa , numOu , numAx
  INTEGER :: i1
  INTEGER :: iIn , iWa , iOu , iAx
  
  numIn = 0 ; numWa = 0 ; numOu = 0 ; numAx = 0
  
  DO i1 =  1 , SIZE(id_Sides)
     IF     ( id_Sides(i1) == 'i' ) THEN
        numIn = numIn + 1
     ELSEIF ( id_Sides(i1) == 'w' ) THEN
        numWa = numWa + 1
     ELSEIF ( id_Sides(i1) == 'o' ) THEN
        numOu = numOu + 1
     ELSEIF ( id_Sides(i1) == 'a' ) THEN
        numAx = numAx + 1
     END IF
  END DO
  ! if some of the numXX is zero, it must be avoided to allocate a vctr of dim = 0
  IF ( numIn .GT. 0) THEN
    ALLOCATE(sides_l_str%inflow (numIn)) ; sides_l_str%inflow  = 0
  ENDIF
  IF ( numWa .GT. 0) THEN
    ALLOCATE(sides_l_str%wall   (numWa)) ; sides_l_str%wall    = 0
  ENDIF
  IF ( numOu .GT. 0) THEN
    ALLOCATE(sides_l_str%outflow(numOu)) ; sides_l_str%outflow = 0
  ENDIF
  IF ( numAx .GT. 0) THEN
    ALLOCATE(sides_l_str%axis   (numAx)) ; sides_l_str%axis    = 0
  ENDIF

  iIn = 0 ; iWa = 0 ; iOu = 0 ; iAx = 0
  DO i1 = 1 , SIZE(id_Sides)
     IF     ( id_Sides(i1) == 'i' .AND. numIn .GT. 0 ) THEN
        iIn = iIn + 1
        sides_l_str%inflow (iIn) = i1
     ELSEIF ( id_Sides(i1) == 'w' .AND. numWa .GT. 0 ) THEN
        iWa = iWa + 1
        sides_l_str%wall   (iWa) = i1
     ELSEIF ( id_Sides(i1) == 'o' .AND. numOu .GT. 0 ) THEN
        iOu = iOu + 1
        sides_l_str%outflow(iOu) = i1
     ELSEIF ( id_Sides(i1) == 'a' .AND. numAx .GT. 0 ) THEN
        iAx = iAx + 1
        sides_l_str%axis   (iAx) = i1
     END IF
  END DO
  
END SUBROUTINE   build_Id_sides_str

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE build_Wal_DirInflow (nCmpn, nos , s_label_str ,  walLogic ,&
                                                            inflLogic )
  
  INTEGER           , INTENT(IN) :: nCmpn , nos
  TYPE(sidesId_str) , INTENT(IN) :: s_label_str
  LOGICAL , DIMENSION(:,:) , ALLOCATABLE , INTENT(INOUT) :: walLogic ,&
                                                            inflLogic
  
  INTEGER :: i1
  
  IF (.NOT. ALLOCATED(walLogic)) THEN
     ALLOCATE(walLogic(nCmpn , nos))
  ENDIF
    IF (.NOT. ALLOCATED(inflLogic)) THEN
     ALLOCATE(inflLogic(nCmpn , nos))
  ENDIF
  
! Build Wal
  walLogic = .FALSE.
  DO i1 = 1 , SIZE(s_label_str%wall)
     walLogic(:,s_label_str%wall(i1)) = .TRUE.
  END DO

! Build DirInflow   
  inflLogic = .FALSE.
  DO i1 = 1 , SIZE(s_label_str%inflow)
     inflLogic(:,s_label_str%inflow(i1)) = .TRUE.
  END DO

   
END SUBROUTINE build_Wal_DirInflow

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE fem_boundary
