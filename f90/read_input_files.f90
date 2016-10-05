MODULE read_input_files

   IMPLICIT NONE

   TYPE program_input

      ! program data
      INTEGER            :: method
      INTEGER            :: nwtn_maxite
      REAL(KIND=8)       :: nwtn_tol
      CHARACTER(LEN=128) :: mesh_directory
      CHARACTER(LEN=128) :: mesh_name
      CHARACTER(LEN=128) :: plot_directory
      CHARACTER(LEN=128) :: restart_directory
      CHARACTER(LEN=128) :: input_restart_file
      CHARACTER(LEN=128) :: output_restart_file
      LOGICAL            :: read_restart_flag
      LOGICAL            :: write_restart_flag
      LOGICAL            :: write_QP_restart_flag
      LOGICAL            :: write_BVS_flag
      LOGICAL            :: write_plots_flag
      ! eigenvalue data
      INTEGER            :: eigen_BC
      INTEGER            :: eigen_nev
      INTEGER            :: eigen_maxit
      REAL(KIND=8)       :: eigen_tol
      COMPLEX(KIND=8)    :: eigen_sigma
      CHARACTER(LEN=128) :: eigen_output_directory
      INTEGER            :: eigen_plotNumber
      INTEGER            :: eigen_directAdjoint_flag
      LOGICAL            :: eigen_compute_structSens_flag
      ! structural sensitivity data
      INTEGER            :: structSens_eigenNumber
      CHARACTER(LEN=128) :: structSens_directEigen_name
      CHARACTER(LEN=128) :: structSens_adjointEigen_name
      ! transient growth data
      INTEGER            :: tranGrowth_method
      INTEGER            :: tranGrowth_initGuess
      INTEGER            :: tranGrowth_BC
      REAL(KIND=8)       :: tranGrowth_tau
      REAL(KIND=8)       :: tranGrowth_dt
      INTEGER            :: tranGrowth_maxit
      REAL(KIND=8)       :: tranGrowth_tol
      ! dns data
      INTEGER            :: dns_method
      REAL(KIND=8)       :: dns_tInit
      REAL(KIND=8)       :: dns_tEnd
      REAL(KIND=8)       :: dns_dt
      REAL(KIND=8)       :: dns_dtPlot
      CHARACTER(LEN=128) :: dns_output_directory
      ! opt perturbation evolution
      CHARACTER(LEN=128) :: etg_opt_perturb
      ! centre manifold reduction
      INTEGER            :: cm_max_order
      REAL(KIND=8)       :: cm_Re
      REAL(KIND=8)       :: cm_tInit
      REAL(KIND=8)       :: cm_tEnd
      REAL(KIND=8)       :: cm_dt 
      ! (9) harmonic gain
      CHARACTER(LEN=128) :: hg_omega_filen
      CHARACTER(LEN=128) :: hg_output_directory
      INTEGER            :: hg_num_compute
      INTEGER            :: hg_num_plot
      CHARACTER(LEN=128) :: hg_forcing_flag
      CHARACTER(LEN=128) :: hg_vol_forcing_domain
      REAL(KIND=8)       :: hg_vol_forcing_blending
      ! (10) stochastic gain
      CHARACTER(LEN=128) :: sg_flag
      ! (11) stochastic gain
      CHARACTER(LEN=128) :: hgSens_flag
      INTEGER            :: hgSens_id
      CHARACTER(LEN=128) :: hgSens_forcing
      CHARACTER(LEN=128) :: hgSens_cyl_flag
      REAL(KIND=8)       :: hgSens_cyl_diam

   END TYPE program_input

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE read_program_data(file_name,  prog_input)

   IMPLICIT NONE
   ! input variables
   CHARACTER(*), INTENT(IN) :: file_name
   ! output variables
   TYPE(program_input) :: prog_input
   ! local variables
   INTEGER :: fid = 22

   ! executable statements
   OPEN (UNIT = fid, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! program data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % method
   READ  (fid,*) prog_input % nwtn_maxite
   READ  (fid,*) prog_input % nwtn_tol
   READ  (fid,*) prog_input % mesh_directory,     prog_input % mesh_name
   READ  (fid,*) prog_input % plot_directory
   READ  (fid,*) prog_input % restart_directory
   READ  (fid,*) prog_input % input_restart_file, prog_input % output_restart_file
   READ  (fid,*) prog_input % read_restart_flag,  prog_input % write_restart_flag
   READ  (fid,*) prog_input % write_QP_restart_flag
   READ  (fid,*) prog_input % write_BVS_flag
   READ  (fid,*) prog_input % write_plots_flag
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (3) eigenvalue data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % eigen_BC
   READ  (fid,*) prog_input % eigen_nev
   READ  (fid,*) prog_input % eigen_maxit
   READ  (fid,*) prog_input % eigen_tol
   READ  (fid,*) prog_input % eigen_sigma
   READ  (fid,*) prog_input % eigen_output_directory
   READ  (fid,*) prog_input % eigen_plotNumber
   READ  (fid,*) prog_input % eigen_directAdjoint_flag
   READ  (fid,*) prog_input % eigen_compute_structSens_flag
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (4) structural sensitivity data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % structSens_eigenNumber
   READ  (fid,*) prog_input % structSens_directEigen_name
   READ  (fid,*) prog_input % structSens_adjointEigen_name
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (5) transient growth data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % tranGrowth_method
   READ  (fid,*) prog_input % tranGrowth_initGuess
   READ  (fid,*) prog_input % tranGrowth_BC
   READ  (fid,*) prog_input % tranGrowth_tau
   READ  (fid,*) prog_input % tranGrowth_dt
   READ  (fid,*) prog_input % tranGrowth_maxit
   READ  (fid,*) prog_input % tranGrowth_tol
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (6) dns data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % dns_method
   READ  (fid,*) prog_input % dns_tInit
   READ  (fid,*) prog_input % dns_tEnd
   READ  (fid,*) prog_input % dns_dt
   READ  (fid,*) prog_input % dns_dtPlot
   READ  (fid,*) prog_input % dns_output_directory
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (7) nonlinear evolution of optimal perturbation
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % etg_opt_perturb
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (8) centre manifold reduction
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % cm_max_order
   READ  (fid,*) prog_input % cm_Re
   READ  (fid,*) prog_input % cm_tInit
   READ  (fid,*) prog_input % cm_tEnd
   READ  (fid,*) prog_input % cm_dt
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (9) harmonic gain
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % hg_omega_filen
   READ  (fid,*) prog_input % hg_output_directory
   READ  (fid,*) prog_input % hg_num_compute
   READ  (fid,*) prog_input % hg_num_plot
   READ  (fid,*) prog_input % hg_forcing_flag
   READ  (fid,*) prog_input % hg_vol_forcing_domain
   READ  (fid,*) prog_input % hg_vol_forcing_blending
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (10) stochastic gain
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % sg_flag
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! (11) stochastic gain
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % hgSens_flag
   READ  (fid,*) prog_input % hgSens_id
   READ  (fid,*) prog_input % hgSens_forcing
   READ  (fid,*) prog_input % hgSens_cyl_flag
   READ  (fid,*) prog_input % hgSens_cyl_diam

   CLOSE(fid)


END SUBROUTINE read_program_data

!==============================================================================

END MODULE read_input_files
