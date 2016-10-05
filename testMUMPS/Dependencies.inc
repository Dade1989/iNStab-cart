###############################################################################
# Objects
#
OBJ    = sparse_matrix_profiles.o    \
         global_variables.o          \
         par_solve_mumps.o
           
###############################################################################
# Dependencies

# Fortran
main.o: $(OBJ)

sparse_matrix_profiles.o:

global_variables.o: sparse_matrix_profiles.o

par_solve_mumps.o: global_variables.o \
                   sparse_matrix_profiles.o

# C
loca_interface.o: loca_wrappers.o
