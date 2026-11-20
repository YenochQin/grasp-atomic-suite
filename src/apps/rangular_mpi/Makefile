EXE=${GRASP}/bin/rangular_mpi
LIBS=-L ${GRASP}/lib/ -lrang90 -lbreit_angular_support -lgrasp_mpi -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/breit_angular_support -I ${GRASP}/src/lib/mpi -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	locate_mcp_restart_position.o \
	locate_mcp_restart_position_interface.o \
	read_calculation_options.o \
	read_calculation_options_interface.o \
	write_mcp_sparsity_summary_mpi.o \
	write_mcp_sparsity_summary_mpi_interface.o \
	configure_debug_output.o \
	configure_debug_output_interface.o \
	configure_debug_output_mpi.o \
	configure_debug_output_mpi_interface.o \
	connect_mcp_files.o \
	connect_mcp_files_interface.o \
	open_mcp_outputs_mpi.o \
	open_mcp_outputs_mpi_interface.o \
	write_mcp_sparsity.o \
	write_mcp_sparsity_interface.o \
	open_summary_file.o \
	open_summary_file_interface.o \
	open_temporary_mcp_channel.o \
	open_temporary_mcp_channel_interface.o \
	sort_mcp_coefficients.o \
	sort_mcp_coefficients_interface.o \
	write_input_summary.o \
	write_input_summary_interface.o \
	compute_mcp_angular_coefficients_mpi.o \
	compute_mcp_angular_coefficients_mpi_interface.o \
	rangular_mpi_main.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
