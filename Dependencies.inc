###############################################################################
# Objects
#
OBJ    = dirichlet_neumann.o         \
         dynamic_structures.o        \
         gauss_points_2d_l.o         \
         gauss_points_2d_p_sp.o      \
         prep_mesh_p1p2_sp.o         \
         qv_sp.o                     \
         qc_sp_m.o                   \
         qs_l_m.o                    \
         qs_l_sp.o                   \
         qs_sp.o                     \
         qs_sp_m.o                   \
         qb_sp_m.o                   \
         sparse_matrix_profiles.o    \
         sparse_matrix_operations.o  \
         global_variables.o          \
         start_sparse_kit.o          \
         fem_2d_grid.o               \
         par_solve_mumps.o           \
         cartesian_boundary_values.o \
         vtk_plot.o                  \
         restart_io.o                \
         read_input_files.o          \
         loca_wrappers.o             \
         loca_interface.o            \
         loca_parameters.o           \
         eigensolve.o                \
         transient_growth.o          \
         dns_algorithms.o            \
         vorticity_stream.o          \
         cm_type_definitions.o       \
         cm_write_complex_arrays.o   \
         cm_sort_module.o            \
         cm_forcing_connectivity.o   \
         cm_index_connectivity.o     \
         cm_forcing_input.o          \
         cm_compute_forcing_terms.o  \
         cm_compute_cm_terms.o       \
         cm_linear.o                 \
         cm_time_evolution.o         \
         centreManifold.o            \
         stressModule.o              \
         hg_build_matrices.o         \
         hg_arpack_module.o          \
         hg_compute_hg_sens.o        \
         harmonicGain.o              \
         fem_boundary.o              \
         fem_start_sparse.o          \
         fem_miscellaneous.o


#fem_start_sparse.o : prep_mesh_p1p2_sp.o       \
#                     global_variables.o        \
#                     sparse_matrix_profiles.o

###############################################################################
# Dependencies

# Fortran
main.o: $(OBJ)

fem_boundary.o : global_variables.o   \
                 fem_miscellaneous.o  \
                 fem_start_sparse.o 

fem_start_sparse.o : prep_mesh_p1p2_sp.o       \
                     global_variables.o        \
                     sparse_matrix_profiles.o  \
                     fem_miscellaneous.o

fem_miscellaneous.o : prep_mesh_p1p2_sp.o  \
                      global_variables.o

hg_compute_hg_sens.o : global_variables.o         \
                       prep_mesh_p1p2_sp.o        \
                       sparse_matrix_profiles.o   \
                       sparse_matrix_operations.o \
                       par_solve_mumps.o          \
                       qc_sp_m.o                  \
                       vtk_plot.o                 \
                       fem_miscellaneous.o

hg_arpack_module.o : global_variables.o          \
                     prep_mesh_p1p2_sp.o         \
                     par_solve_mumps.o           \
                     sparse_matrix_operations.o  \
                     cartesian_boundary_values.o \
                     vtk_plot.o                  \
                     fem_miscellaneous.o

hg_build_matrices.o : par_solve_mumps.o           \
                      sparse_matrix_profiles.o    \
                      sparse_matrix_operations.o  \
                      start_sparse_kit.o          \
                      fem_boundary.o              \
                      fem_start_sparse.o          \
                      fem_miscellaneous.o         \
                      fem_start_sparse.o          \
                      qb_sp_m.o                   \
                      prep_mesh_p1p2_sp.o         \
                      global_variables.o

harmonicGain.o : hg_build_matrices.o         \
                 hg_arpack_module.o          \
                 hg_compute_hg_sens.o        \
                 sparse_matrix_profiles.o    \
                 global_variables.o

stressModule.o : sparse_matrix_profiles.o     \
                 sparse_matrix_operations.o   \
                 global_variables.o           \
                 prep_mesh_p1p2_sp.o          \
                 gauss_points_2d_p_sp.o       \
                 gauss_points_2d_l.o          \
                 qb_sp_m.o                    \
                 par_solve_mumps.o

cm_time_evolution.o : global_variables.o

cm_type_definitions.o :

cm_write_complex_arrays.o :

cm_sort_module.o :

cm_forcing_connectivity.o : cm_type_definitions.o

cm_index_connectivity : cm_type_definitions.o \
                        cm_sort_module.o

cm_forcing_input.o : sparse_matrix_operations.o    \
                     sparse_matrix_profiles.o      \
                     global_variables.o            \
                     start_sparse_kit.o            \
                     gauss_points_2d_p_sp.o        \
                     prep_mesh_p1p2_sp.o           \
                     qc_sp_m.o                     \
                     qv_sp.o                       \
                     vtk_plot.o

cm_compute_cm_terms.o : cm_compute_forcing_terms.o   \
                        cm_index_connectivity.o      \
                        cm_linear.o                  \
                        cm_write_complex_arrays.o    \
                        global_variables.o           \
                        eigensolve.o                 \
                        prep_mesh_p1p2_sp.o          \
                        stressModule.o

cm_compute_forcing_terms.o : global_variables.o        \
                             cm_forcing_connectivity.o \
                             cm_index_connectivity.o   \
                             cm_type_definitions.o     \
                             cm_forcing_input.o

cm_linear.o : global_variables.o        \
              sparse_matrix_profiles.o  \
              dirichlet_neumann.o       \
              cartesian_boundary_values.o \
              qc_sp_m.o                 \
              eigensolve.o              \
              restart_io.o

centreManifold.o: global_variables.o        \
                  eigensolve.o              \
                  restart_io.o              \
                  sparse_matrix_profiles.o  \
                  par_solve_mumps.o         \
                  cm_forcing_input.o        \
                  cm_linear.o               \
                  cm_index_connectivity.o   \
                  cm_compute_cm_terms.o

qb_sp_m.o: sparse_matrix_profiles.o \
           global_variables.o
 
read_input_files.o:
dirichlet_neumann.o:
dynamic_structures.o:
gauss_points_2d_l.o:
gauss_points_2d_p_sp.o:
fem_2d_grid.o:
sparse_matrix_profiles.o:
sparse_matrix_operations.o: global_variables.o   \
                            dirichlet_neumann.o  \
                            prep_mesh_p1p2_sp.o
loca_parameters.o:
global_variables.o: read_input_files.o         \
                    dynamic_structures.o       \
                    sparse_matrix_profiles.o

par_solve_mumps.o: global_variables.o \
                   sparse_matrix_profiles.o

start_sparse_kit.o: sparse_matrix_profiles.o \
                    prep_mesh_p1p2_sp.o

qs_l_m.o: sparse_matrix_profiles.o \
          gauss_points_2d_p_sp.o   \
          gauss_points_2d_l.o

qs_l_sp.o: gauss_points_2d_p_sp.o

qs_sp.o: gauss_points_2d_p_sp.o

qs_sp_m.o: sparse_matrix_profiles.o \
           gauss_points_2d_p_sp.o

qv_sp.o: gauss_points_2d_p_sp.o

qc_sp_m.o: sparse_matrix_profiles.o \
           gauss_points_2d_l.o      \
           gauss_points_2d_p_sp.o

prep_mesh_p1p2_sp.o: gauss_points_2d_l.o    \
                     gauss_points_2d_p_sp.o \
                     dirichlet_neumann.o    \
                     fem_2d_grid.o

vtk_plot.o: start_sparse_kit.o

restart_io.o: sparse_matrix_profiles.o \
              global_variables.o       \
              prep_mesh_p1p2_sp.o      \
              start_sparse_kit.o

cartesian_boundary_values.o: global_variables.o   \
                             dynamic_structures.o \
                             dirichlet_neumann.o  \
                             gauss_points_2d_p_sp.o

loca_wrappers.o: loca_parameters.o           \
                 par_solve_mumps.o           \
                 qv_sp.o                     \
                 qc_sp_m.o                   \
                 qs_l_sp.o                   \
                 dynamic_structures.o        \
                 sparse_matrix_profiles.o    \
                 global_variables.o          \
                 prep_mesh_p1p2_sp.o         \
                 start_sparse_kit.o          \
                 dirichlet_neumann.o         \
                 read_input_files.o          \
                 eigensolve.o                \
                 cartesian_boundary_values.o \
                 vtk_plot.o

eigensolve.o: sparse_matrix_profiles.o   \
              sparse_matrix_operations.o \
              start_sparse_kit.o

transient_growth.o: dynamic_structures.o       \
                    sparse_matrix_profiles.o   \
                    sparse_matrix_operations.o \
                    global_variables.o         \
                    prep_mesh_p1p2_sp.o        \
                    start_sparse_kit.o         \
                    dirichlet_neumann.o        \
                    qc_sp_m.o                  \
                    qv_sp.o                    \
                    qs_l_m.o                   \
                    qs_l_sp.o                  \
                    par_solve_mumps.o          \
                    vtk_plot.o

dns_algorithms.o: dynamic_structures.o     \
                  sparse_matrix_profiles.o \
                  global_variables.o       \
                  prep_mesh_p1p2_sp.o      \
                  start_sparse_kit.o       \
                  qc_sp_m.o                \
                  qv_sp.o                  \
                  par_solve_mumps.o        \
                  restart_io.o             \
                  vtk_plot.o

vorticity_stream.o: global_variables.o  \
                    prep_mesh_p1p2_sp.o \
                    dirichlet_neumann.o \
                    start_sparse_kit.o  \
                    qs_sp.o             \
                    qs_sp_m.o           \
                    par_solve_mumps.o

# C
loca_interface.o: loca_wrappers.o
