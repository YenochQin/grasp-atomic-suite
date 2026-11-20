EXE=rangular
LIBRARIES="spin_angular breit_angular_support core grasp_shared"
LAPACK=true
FILES="
try_allocate_sort_buffers.f90
locate_mcp_restart_position.f90
read_calculation_options.f90
sort_mcp_coefficients.f90
sort_mcp_coefficients_in_memory.f90
write_mcp_sparsity_summary.f90
write_mcp_sparsity.f90
compute_mcp_angular_coefficients.f90
configure_debug_output.f90
connect_mcp_files.f90
open_mcp_outputs.f90
open_summary_file.f90
open_temporary_mcp_channels.f90
write_input_summary.f90

rangular_main.f90

# cons_C was not being compiled in the original makefile
#numeric_constants.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
