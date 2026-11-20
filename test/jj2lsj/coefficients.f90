! Exercise the real parser, LS basis generator and transformation, then expose
! every signed double for bitwise comparison, including the zero truncated tail.
program coefficients
    use jj2lsj_code, only: asf2ls, asf_set_ls, blocks_number, dallocasfls, epsnew, inscreen, setls, vectors_number
    use orbital_csf_basis, only: NCF => csf_count
    use selected_level_indices, only: NVEC
    use eigenstate_expansions, only: EVEC => csf_mixing_coefficients
    use transition_basis_blocks, only: NBLOCK, NCFINBLK
    use log_factorial_table, only: factt
    use physical_constants, only: setcon
    use isotope_data_input, only: setiso
    use, intrinsic :: iso_fortran_env, only: real64, int64
    implicit none (type, external)
    integer :: coefficient_output_unit
    integer :: counts(Blocks_number)
    integer, allocatable :: levels(:,:)
    integer :: block, first, last, lev, level, i, outc, outj, unique, outt
    integer, allocatable :: mask(:)
    integer(int64) :: stride
    real(real64) :: threshold, weight

    call setcon
    call factt
    call setiso("isodata")
    allocate(levels(Blocks_number,Vectors_number))
    call inscreen(threshold,levels,counts,outc,outj,unique,outt)
    allocate(mask(NCF))
    stride = NCF
    last = 0
    open(newunit=coefficient_output_unit,file="coefficients.raw",status="replace",access="stream",form="unformatted")
    do block = 1,NBLOCK
        first = last+1
        last = first+NCFINBLK(block)-1
        if (counts(block) == 0) cycle
        mask = 0
        do i = first,last
            weight = 0.0_real64
            do lev = 1,counts(block)
                level = levels(block,lev)
                weight = weight+abs(EVEC(i+(level-1)*stride))
            end do
            if (abs(weight) >= abs(EPSNEW)) mask(i) = 1
        end do
        call setLS(mask,first,last)
        allocate(asf_set_LS%asf(NVEC))
        do i = 1,NVEC
            allocate(asf_set_LS%asf(i)%eigenvector(asf_set_LS%csf_set_LS%nocsf))
            asf_set_LS%asf(i)%eigenvector = 0.0_real64
        end do
        call asf2ls(1,mask,counts(block),block,levels,first,last,outt)
        write(coefficient_output_unit) block, counts(block), asf_set_LS%csf_set_LS%nocsf
        do lev = 1,counts(block)
            write(coefficient_output_unit) asf_set_LS%asf(levels(block,lev))%eigenvector
        end do
        call dallocASFLS(asf_set_LS)
    end do
    close(coefficient_output_unit)
    print *, "Coefficient dump complete"
end program coefficients
