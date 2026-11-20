EXE=${GRASP}/bin/rangular
LIBS=-L ${GRASP}/lib/ -lrang90 -lbreit_angular_support -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/breit_angular_support -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	try_allocate_sort_buffers.o \
	try_allocate_sort_buffers_interface.o \
	locate_mcp_restart_position.o \
	locate_mcp_restart_position_interface.o \
	read_calculation_options.o \
	read_calculation_options_interface.o \
	sort_mcp_coefficients.o \
	sort_mcp_coefficients_interface.o \
	sort_mcp_coefficients_in_memory.o \
	sort_mcp_coefficients_in_memory_interface.o \
	write_mcp_sparsity_summary.o \
	write_mcp_sparsity_summary_interface.o \
	write_mcp_sparsity.o \
	write_mcp_sparsity_interface.o \
	compute_mcp_angular_coefficients.o \
	compute_mcp_angular_coefficients_interface.o \
	configure_debug_output.o \
	configure_debug_output_interface.o \
	connect_mcp_files.o \
	connect_mcp_files_interface.o \
	open_mcp_outputs.o \
	open_mcp_outputs_interface.o \
	open_summary_file.o \
	open_summary_file_interface.o \
	open_temporary_mcp_channels.o \
	open_temporary_mcp_channels_interface.o \
	write_input_summary.o \
	write_input_summary_interface.o \
	rangular_main.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
