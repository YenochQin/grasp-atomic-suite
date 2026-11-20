EXE=rangular_mpi
LIBRARIES="spin_angular breit_angular_support grasp_mpi core grasp_shared"
LAPACK=true
ISMPI=true
FILES="
locate_mcp_restart_position.f90
read_calculation_options.f90
write_mcp_sparsity_summary_mpi.f90
configure_debug_output.f90
configure_debug_output_mpi.f90
connect_mcp_files.f90
open_mcp_outputs_mpi.f90
write_mcp_sparsity.f90
open_summary_file.f90
open_temporary_mcp_channel.f90
sort_mcp_coefficients.f90
write_input_summary.f90
compute_mcp_angular_coefficients_mpi.f90
rangular_mpi_main.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
