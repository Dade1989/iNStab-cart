MODULE global_variables

   PUBLIC

! global variables for MUMPS 

   INCLUDE 'mpif.h'
   INTEGER :: myRank, nProc, mpiIerr
   INTEGER :: parMumpsJob
   INTEGER :: matrID

END MODULE global_variables
