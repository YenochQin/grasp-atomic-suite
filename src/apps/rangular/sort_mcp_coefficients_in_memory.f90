module sort_mcp_coefficients_in_memory_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: sort_mcp_coefficients_in_memory
contains

    !***********************************************************************
    !                                                                      *
    subroutine sort_mcp_coefficients_in_memory(mcp_channel_index,sort_capacity_per_channel,coefficient_count, &
        & distinct_integral_count,print_details, &
        & symmetry_block_index,output_file_prefix)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   This routine sorts lists                                           *
        !                                                                      *
        !                     (ICLMN,INDEX,LABEL,COEFF)                        *
        !   into lists                                                         *
        !                     (LABEL,ICLMN,INDEX,COEFF)                        *
        !                                                                      *
        !   using Heapsort. File NFILE is closed by this routine.  NCOEFF is   *
        !   the number of triads (INDEX, ...) and is an input. NTGRAL is the   *
        !   number of different values of LABEL, and is an output. If LPRINT   *
        !   is  .TRUE. , the contents of the sorted file are interpreted and   *
        !   printed to unit 99.                                                *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use array_memory_management, only: alloc, dalloc
        use orbital_csf_basis,           only: orbital_principal_n, csf_count, orbital_symmetry_labels
        use standard_io_units, only: stderr_unit
        use angular_coefficient_channels, only: channel_coefficients, channel_packed_orbital_labels, &
            & channel_csf_columns, channel_matrix_positions
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Local declarations
        integer :: mcp_append_unit
        integer :: heap_parent_index, heap_child_index

        ! Arguments
        integer, intent(in) :: mcp_channel_index
        integer, intent(in) :: coefficient_count
        integer, intent(out)  :: distinct_integral_count
        integer, intent(in)   :: symmetry_block_index,sort_capacity_per_channel
        logical, intent(in)   :: print_details
        character, intent(in) :: output_file_prefix*(*)

        ! Local constants

        ! Local declarations
        !      INTEGER, DIMENSION(:), pointer :: ICLMN, INDEX, LABEL, NSWAP
        integer, pointer :: nswap(:)
        !      REAL(DOUBLE), DIMENSION(:), pointer :: COEFF
        integer :: text_length, lck, io_status, i, heap_root_index, heap_end_index, icl, ind, lab, nsw, j, last, ibeg, iend, &
            & ncontr, ia, ib, k, id, ic, &
            igg
        real(real64) :: cof
        character :: number_text*20, ck*2
        character (len = len (output_file_prefix) + 3):: file_path
        !-----------------------------------------------
        !
        if(mcp_channel_index == 31) then
            igg = 1
        else
            igg = mcp_channel_index - 30
        end if
        !      IGG = IGGMAX_K*(IGG-1)
        !   Message
        !
        call format_integer (coefficient_count, number_text, text_length)
        if (mcp_channel_index > 31) then
            call format_integer (mcp_channel_index - 32, ck, lck)
            write (output_unit, *) "Sorting "//number_text(1:text_length)//" V(k="//ck(1:lck)//&
                ") coefficients ...", mcp_channel_index
        else
            write (output_unit, *) "Sorting "//number_text(1:text_length)//" T coefficients ...", mcp_channel_index
        end if
        !
        call format_integer (mcp_channel_index, ck, lck)
        if (lck > 2) then
            write (stderr_unit, *) "sort: nfile > 99; check fullname"
            stop
        end if

        file_path = output_file_prefix//"."//ck(1:2)

        open(newunit=mcp_append_unit,file=file_path,status="OLD",form="UNFORMATTED", &
            iostat=io_status, position="APPEND")
        if (io_status /= 0) then
            write (stderr_unit, *) " Error when opening the file ", file_path
            stop
        end if
        !
        !   Sort the list
        !
        if (coefficient_count > 0) then
            !
            !   Allocate storage for all required arrays
            !
            if (mcp_channel_index == 33) call alloc (nswap, coefficient_count,"NSWAP", "SORTMEM")
        end if
        !
        !   Sort LABEL into ascending order using the heapsort algorithm;
        !   move the associated members of COEFF and INDEX in the same
        !   manner; the code below is adapted from Press et al.
        !
        if (mcp_channel_index==33 .and. coefficient_count>1) then

            heap_root_index = coefficient_count/2 + 1
            heap_end_index = coefficient_count
            do
                if (heap_root_index > 1) then
                    heap_root_index = heap_root_index - 1
                    cof = channel_coefficients(heap_root_index,igg)
                    icl = channel_csf_columns(heap_root_index,igg)
                    ind = channel_matrix_positions(heap_root_index,igg)
                    lab = channel_packed_orbital_labels(heap_root_index,igg)
                    nsw = nswap(heap_root_index)
                else
                    cof = channel_coefficients(heap_end_index,igg)
                    icl = channel_csf_columns(heap_end_index,igg)
                    ind = channel_matrix_positions(heap_end_index,igg)
                    lab = channel_packed_orbital_labels(heap_end_index,igg)
                    nsw = nswap(heap_end_index)
                    channel_coefficients(heap_end_index,igg) = channel_coefficients(1,igg)
                    channel_csf_columns(heap_end_index,igg) = channel_csf_columns(1,igg)
                    channel_matrix_positions(heap_end_index,igg) = channel_matrix_positions(1,igg)
                    channel_packed_orbital_labels(heap_end_index,igg) = channel_packed_orbital_labels(1,igg)
                    nswap(heap_end_index) = nswap(1)
                    heap_end_index = heap_end_index - 1
                    if (heap_end_index == 1) then
                        channel_coefficients(1,igg) = cof
                        channel_csf_columns(1,igg) = icl
                        channel_matrix_positions(1,igg) = ind
                        channel_packed_orbital_labels(1,igg) = lab
                        nswap(1) = nsw
                        exit
                    end if
                end if
                heap_parent_index = heap_root_index
                heap_child_index = heap_root_index + heap_root_index
                do while (heap_child_index <= heap_end_index)
                    if (heap_child_index < heap_end_index) then
                        if (channel_packed_orbital_labels(heap_child_index,igg) < &
                            & channel_packed_orbital_labels(heap_child_index+1,igg)) heap_child_index = heap_child_index + 1
                    end if
                    if (lab < channel_packed_orbital_labels(heap_child_index,igg)) then
                        channel_coefficients(heap_parent_index,igg) = channel_coefficients(heap_child_index,igg)
                        channel_csf_columns(heap_parent_index,igg) = channel_csf_columns(heap_child_index,igg)
                        channel_matrix_positions(heap_parent_index,igg) = channel_matrix_positions(heap_child_index,igg)
                        channel_packed_orbital_labels(heap_parent_index,igg) = channel_packed_orbital_labels(heap_child_index,igg)
                        nswap(heap_parent_index) = nswap(heap_child_index)
                        heap_parent_index = heap_child_index
                        heap_child_index = heap_child_index + heap_child_index
                    else
                        heap_child_index = heap_end_index + 1
                    end if
                end do
                channel_coefficients(heap_parent_index,igg) = cof
                channel_csf_columns(heap_parent_index,igg) = icl
                channel_matrix_positions(heap_parent_index,igg) = ind
                channel_packed_orbital_labels(heap_parent_index,igg) = lab
                nswap(heap_parent_index) = nsw
            end do

        else if (mcp_channel_index/=33 .and. coefficient_count>1) then
            !
            !   Sort LABEL into ascending order using the heapsort algorithm;
            !   move the associated members of COEFF and INDEX in the same
            !   manner; the code below is adapted from Press et al.
            !
            heap_root_index = coefficient_count/2 + 1
            heap_end_index = coefficient_count
            do
                if (heap_root_index > 1) then
                    heap_root_index = heap_root_index - 1
                    cof = channel_coefficients(heap_root_index,igg)
                    icl = channel_csf_columns(heap_root_index,igg)
                    ind = channel_matrix_positions(heap_root_index,igg)
                    lab = channel_packed_orbital_labels(heap_root_index,igg)
                else
                    cof = channel_coefficients(heap_end_index,igg)
                    icl = channel_csf_columns(heap_end_index,igg)
                    ind = channel_matrix_positions(heap_end_index,igg)
                    lab = channel_packed_orbital_labels(heap_end_index,igg)
                    channel_coefficients(heap_end_index,igg) = channel_coefficients(1,igg)
                    channel_csf_columns(heap_end_index,igg) = channel_csf_columns(1,igg)
                    channel_matrix_positions(heap_end_index,igg) = channel_matrix_positions(1,igg)
                    channel_packed_orbital_labels(heap_end_index,igg) = channel_packed_orbital_labels(1,igg)
                    heap_end_index = heap_end_index - 1
                    if (heap_end_index == 1) then
                        channel_coefficients(1,igg) = cof
                        channel_csf_columns(1,igg) = icl
                        channel_matrix_positions(1,igg) = ind
                        channel_packed_orbital_labels(1,igg) = lab
                        exit
                    end if
                end if
                heap_parent_index = heap_root_index
                heap_child_index = heap_root_index + heap_root_index
                do while (heap_child_index <= heap_end_index)
                    if (heap_child_index < heap_end_index) then
                        if (channel_packed_orbital_labels(heap_child_index,igg) < &
                            & channel_packed_orbital_labels(heap_child_index+1,igg)) heap_child_index = heap_child_index + 1
                    end if
                    if (lab < channel_packed_orbital_labels(heap_child_index,igg)) then
                        channel_coefficients(heap_parent_index,igg) = channel_coefficients(heap_child_index,igg)
                        channel_csf_columns(heap_parent_index,igg) = channel_csf_columns(heap_child_index,igg)
                        channel_matrix_positions(heap_parent_index,igg) = channel_matrix_positions(heap_child_index,igg)
                        channel_packed_orbital_labels(heap_parent_index,igg) = channel_packed_orbital_labels(heap_child_index,igg)
                        heap_parent_index = heap_child_index
                        heap_child_index = heap_child_index + heap_child_index
                    else
                        heap_child_index = heap_end_index + 1
                    end if
                end do
                channel_coefficients(heap_parent_index,igg) = cof
                channel_csf_columns(heap_parent_index,igg) = icl
                channel_matrix_positions(heap_parent_index,igg) = ind
                channel_packed_orbital_labels(heap_parent_index,igg) = lab
            end do

        end if
        !
        !   Sorting complete; rewrite the file header
        !
        !
        write (mcp_append_unit) "MCP", symmetry_block_index, csf_count, coefficient_count
        !GG      WRITE (9999,*) 'MCP', NB, NCF, NCOEFF, NFILE
        !
        !   Write the sorted list to mcp.xx
        !
        if (coefficient_count > 0) then
            !
            last = channel_packed_orbital_labels(1,igg)
            ibeg = 1
            iend = 1
            distinct_integral_count = 1
            !
            do i = 2, coefficient_count
                if (channel_packed_orbital_labels(i,igg) == last) then
                    iend = iend + 1
                else
                    write (mcp_append_unit) last, iend - ibeg + 1
                    !GG               WRITE (9999,*) LAST, IEND - IBEG + 1
                    write (mcp_append_unit) (channel_csf_columns(j,igg),channel_matrix_positions(j,igg),channel_coefficients( &
                        & j,igg), &
                        j=ibeg,iend)
                    !cjb-GG              DO J = IBEG,IEND
                    !GG               WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
                    !cjb-GG              END DO

                    !              IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,IEND)
                    distinct_integral_count = distinct_integral_count + 1
                    last = channel_packed_orbital_labels(i,igg)
                    ibeg = iend + 1
                    iend = ibeg
                end if
            end do
            !
            if (ibeg <= coefficient_count) then
                write (mcp_append_unit) last, coefficient_count - ibeg + 1
                !GG            WRITE (9999,*) LAST, NCOEFF - IBEG + 1
                write (mcp_append_unit) (channel_csf_columns(j,igg),channel_matrix_positions(j,igg),channel_coefficients(j, &
                    & igg),       &
                    j=ibeg,coefficient_count)
                !GG          DO J=IBEG,NCOEFF
                !GG            WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
                !GG          END DO

                !           IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,NCOEFF)
            end if
        !
        else
            !
            distinct_integral_count = 0
            !
        end if
        !
        !     write the terminator record for this block
        !
        write (mcp_append_unit) 0, 0
        !GG      WRITE (9999,*) 0, 0
        close(mcp_append_unit)
        !
        !   Completion message
        !
        !PRINT *, ' ... sort complete; ', ntgral, ' integrals;'
        !
        !   Debug printout
        !
        if (print_details) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"From SORT:")')
            write (output_unit, '(/,"From SORT:")')
            if (coefficient_count > 0) then
                !
                last = channel_packed_orbital_labels(1,igg)
                ibeg = 1
                iend = 1
                !
                do i = 2, coefficient_count
                    if (channel_packed_orbital_labels(i,igg) == last) iend = iend + 1
                    do
                        if (channel_packed_orbital_labels(i,igg)==last .and. i/=coefficient_count) exit
                        lab = last
                        ncontr = iend - ibeg + 1
                        if (mcp_channel_index == 31) then
                            ia = mod(lab,orbital_index_radix)
                            ib = lab/orbital_index_radix
                            call ensure_debug_output()
                            write (angular_debug_unit, '(" I(",1i2,1a2,",",1i2,1a2,"):")') &
                                orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                                & orbital_symmetry_labels(ib)
                            do j = ibeg, iend
                                call ensure_debug_output()
                                write (angular_debug_unit, '("  T_[",1i2,",",1i4,"] = ",1p,d19.12)') channel_csf_columns(j,igg), &
                                    & channel_matrix_positions(j,igg),       &
                                    channel_coefficients(j,igg)
                            end do
                        else
                            k = mcp_channel_index - 32
                            id = mod(lab,orbital_index_radix)
                            lab = lab/orbital_index_radix
                            ib = mod(lab,orbital_index_radix)
                            lab = lab/orbital_index_radix
                            ic = mod(lab,orbital_index_radix)
                            ia = lab/orbital_index_radix
                            call ensure_debug_output()
                            write (angular_debug_unit, '(" R^[(",1i2,")] (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2,1a2,"):")') &
                                k, orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                                & orbital_symmetry_labels(ib), orbital_principal_n(ic), orbital_symmetry_labels&
                                (ic), orbital_principal_n(id), orbital_symmetry_labels(id)
                            do j = ibeg, iend
                                call ensure_debug_output()
                                write (angular_debug_unit, '("  V^[(",1i2,")]_[",1i8,",",1i8,"] = ",1p,d19.12)') &
                                    k,channel_csf_columns(j,igg),channel_matrix_positions(j,igg),     &
                                    channel_coefficients(j,igg)
                            end do
                        end if
                        last = channel_packed_orbital_labels(i,igg)
                        ibeg = iend + 1
                        iend = ibeg
                        if (iend /= coefficient_count) exit
                    end do
                end do
            end if
            call ensure_debug_output()
            write (angular_debug_unit, '("  Number of integrals is ",1i4)') distinct_integral_count
        end if
        !
        !   Deallocate storage
        !
        if (coefficient_count > 0) then
            if (mcp_channel_index == 33) call dalloc (nswap, "NSWAP", "SORTMEM")
        end if


        return
    end subroutine sort_mcp_coefficients_in_memory
end module sort_mcp_coefficients_in_memory_interface
