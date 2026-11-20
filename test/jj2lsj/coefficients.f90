! Exercise the real parser, LS basis generator and transformation, then expose
! every signed double for bitwise comparison, including the zero truncated tail.
program coefficients
    use jj2lsj_conversion, only: transform_asf_to_ls, asf_set_ls, max_symmetry_blocks, release_ls_expansion, &
        & coefficient_threshold, &
       & read_ls_conversion_options, build_ls_basis, max_selected_levels
    use orbital_csf_basis, only: csf_count
    use selected_level_indices, only: active_level_count
    use eigenstate_expansions, only: csf_mixing_coefficients
    use transition_basis_blocks, only: transition_block_count, transition_csf_counts_by_block
    use log_factorial_table, only: initialize_log_factorial_table
    use physical_constants, only: initialize_physical_constants
    use isotope_data_input, only: read_isotope_data
    use, intrinsic :: iso_fortran_env, only: real64, int64
    implicit none (type, external)
    integer :: coefficient_output_unit
    integer :: counts(max_symmetry_blocks)
    integer, allocatable :: levels(:,:)
    integer :: block, first, last, lev, level, i, outc, outj, unique, outt
    integer, allocatable :: mask(:)
    integer(int64) :: stride
    real(real64) :: threshold, weight

    call initialize_physical_constants
    call initialize_log_factorial_table
    call read_isotope_data("isodata")
    allocate(levels(max_symmetry_blocks,max_selected_levels))
    call read_ls_conversion_options(threshold,levels,counts,outc,outj,unique,outt)
    allocate(mask(csf_count))
    stride = csf_count
    last = 0
    open(newunit=coefficient_output_unit,file="coefficients.raw",status="replace",access="stream",form="unformatted")
    do block = 1,transition_block_count
        first = last+1
        last = first+transition_csf_counts_by_block(block)-1
        if (counts(block) == 0) cycle
        mask = 0
        do i = first,last
            weight = 0.0_real64
            do lev = 1,counts(block)
                level = levels(block,lev)
                weight = weight+abs(csf_mixing_coefficients(i+(level-1)*stride))
            end do
            if (abs(weight) >= abs(coefficient_threshold)) mask(i) = 1
        end do
        call build_ls_basis(mask,first,last)
        allocate(asf_set_LS%asf(active_level_count))
        do i = 1,active_level_count
            allocate(asf_set_LS%asf(i)%eigenvector(asf_set_LS%csf_set_LS%ls_csf_count))
            asf_set_LS%asf(i)%eigenvector = 0.0_real64
        end do
        call transform_asf_to_ls(1,mask,counts(block),block,levels,first,last,outt)
        write(coefficient_output_unit) block, counts(block), asf_set_LS%csf_set_LS%ls_csf_count
        do lev = 1,counts(block)
            write(coefficient_output_unit) asf_set_LS%asf(levels(block,lev))%eigenvector
        end do
        call release_ls_expansion(asf_set_LS)
    end do
    close(coefficient_output_unit)
    print *, "Coefficient dump complete"
end program coefficients
