!> \brief crea vettori e matrici di lunghezza variabile
!!        semplicemente definendo i nuovi tipi scritti qui sotto
!!        usando l'attributo POINTER invece di ALLOCATABLE
!! 
!!  MODULE dynamic_structures
!!
!!      TYPE dyn_int_line
!!          INTEGER,      DIMENSION(:),   POINTER :: DIL
!!      END TYPE dyn_int_line
!!
!!      TYPE dyn_real_line
!!          REAL(KIND=8), DIMENSION(:),   POINTER :: DRL
!!      END TYPE dyn_real_line
!!
!!      TYPE dyn_real_array
!!          REAL(KIND=8), DIMENSION(:,:), POINTER :: DRA
!!      END TYPE dyn_real_array
!!
!!  END MODULE dynamic_structures
MODULE dynamic_structures



   TYPE dyn_int_line
     
      INTEGER, DIMENSION(:), POINTER :: DIL  ! POINTER attribute
   
   END TYPE dyn_int_line                     ! instead of ALLOCATABLE



   TYPE dyn_real_line

      REAL(KIND=8), DIMENSION(:), POINTER :: DRL ! POINTER attribute

   END TYPE dyn_real_line                        ! instead of ALLOCATABLE


   TYPE dyn_real_array

      REAL(KIND=8), DIMENSION(:,:), POINTER :: DRA ! POINTER attribute

   END TYPE dyn_real_array                         ! instead of ALLOCATABLE



END MODULE dynamic_structures
