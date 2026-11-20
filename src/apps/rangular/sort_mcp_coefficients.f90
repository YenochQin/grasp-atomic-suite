!***********************************************************************
!                                                                      *
subroutine sort_mcp_coefficients(mcp_channel_index, coefficient_count, distinct_integral_count, print_details, &
    & symmetry_block_index, &
    & output_file_prefix)
    ! Module dependencies
    use mcp_generation_options, only: mcp_temporary_units
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
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use integer_format, only: convrt
    implicit none (type, external)

    ! Local declarations
    integer :: mcp_append_unit

    ! Arguments
    integer, intent(in) :: mcp_channel_index
    integer, intent(in) :: coefficient_count
    integer, intent(out)  :: distinct_integral_count
    integer, intent(in)   :: symmetry_block_index
    logical, intent(in)   :: print_details
    character, intent(in) :: output_file_prefix*(*)

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix

    ! Local declarations
    integer, pointer :: iclmn(:), index(:), label(:), nswap(:)
    real(real64), pointer :: coeff(:)
    integer :: text_length, lck, io_status, mb, i, l, ir, icl, ind, lab, nsw, j,&
        last, ibeg, iend, ncontr, ia, ib, k, id, ic
    real(real64) :: cof
    character :: number_text*20, srtlab*8, mcp_channel_label*3, ck*2
    character (len = len (output_file_prefix) + 3):: file_path
    !-----------------------------------------------
    !
    !
    !     POINTER (PCOEFF,COEFF(1)), (PICLMN,ICLMN(1)), (PINDEX,INDEX(1)),
    !    : (PLABEL,LABEL(1)), (PNSWAP,NSWAP(1))
    !
    !     POINTER (PNTRIQ,RIQDUMMY)
    !
    !
    !   This is the encoding key
    !
    !   Message
    !
    call convrt (coefficient_count, number_text, text_length)
    if (mcp_channel_index > 31) then
        call convrt (mcp_channel_index - 32, ck, lck)
        write (output_unit, *) "Sorting "//number_text(1:text_length)//" V(k="//ck(1:lck)//&
            ") coefficients ...", mcp_channel_index
    else
        write (output_unit, *) "Sorting "//number_text(1:text_length)//" T coefficients ...", mcp_channel_index
    end if
    !
    call convrt (mcp_channel_index, ck, lck)
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

    !     IF (NFILE.EQ.33) REWIND (20)

    read (mcp_temporary_units(mcp_channel_index)) mcp_channel_label, mb
    if (symmetry_block_index /= mb) then
        write (stderr_unit, *) "sort: nb = ", symmetry_block_index, ".NE. mb (=", mb, ")"
        stop
    end if
    !
    !   Sort the list
    !
    if (coefficient_count > 0) then
        !
        !   Allocate storage for all required arrays
        !
        call alloc (coeff, coefficient_count, "COEFF", "SORT")
        call alloc (iclmn, coefficient_count, "ICLMN", "SORT")
        call alloc (index, coefficient_count, "INDEX", "SORT")
        call alloc (label, coefficient_count, "LABEL", "SORT")

        if (mcp_channel_index == 33) call alloc (nswap, coefficient_count,"NSWAP", "SORT"  )
        !
        !   Read arrays into memory from NFILE
        !
        do i = 1, coefficient_count
            read (mcp_temporary_units(mcp_channel_index)) iclmn(i), index(i), label(i), coeff(i)
        end do

        !         IF (NFILE.EQ.33) THEN
        !            DO 11 I = 1,NCOEFF
        !               READ (20) NSWAP(I)
        !   11       CONTINUE
        !         ENDIF
        !
    end if
    !
    !   Sort LABEL into ascending order using the heapsort algorithm;
    !   move the associated members of COEFF and INDEX in the same
    !   manner; the code below is adapted from Press et al.
    !
    if (mcp_channel_index==33 .and. coefficient_count>1) then

        l = coefficient_count/2 + 1
        ir = coefficient_count
        do
            if (l > 1) then
                l = l - 1
                cof = coeff(l)
                icl = iclmn(l)
                ind = index(l)
                lab = label(l)
                nsw = nswap(l)
            else
                cof = coeff(ir)
                icl = iclmn(ir)
                ind = index(ir)
                lab = label(ir)
                nsw = nswap(ir)
                coeff(ir) = coeff(1)
                iclmn(ir) = iclmn(1)
                index(ir) = index(1)
                label(ir) = label(1)
                nswap(ir) = nswap(1)
                ir = ir - 1
                if (ir == 1) then
                    coeff(1) = cof
                    iclmn(1) = icl
                    index(1) = ind
                    label(1) = lab
                    nswap(1) = nsw
                    exit
                end if
            end if
            i = l
            j = l + l
            do while (j <= ir)
                if (j < ir) then
                    if (label(j) < label(j+1)) j = j + 1
                end if
                if (lab < label(j)) then
                    coeff(i) = coeff(j)
                    iclmn(i) = iclmn(j)
                    index(i) = index(j)
                    label(i) = label(j)
                    nswap(i) = nswap(j)
                    i = j
                    j = j + j
                else
                    j = ir + 1
                end if
            end do
            coeff(i) = cof
            iclmn(i) = icl
            index(i) = ind
            label(i) = lab
            nswap(i) = nsw
        end do

    else if (mcp_channel_index/=33 .and. coefficient_count>1) then
        !
        !   Sort LABEL into ascending order using the heapsort algorithm;
        !   move the associated members of COEFF and INDEX in the same
        !   manner; the code below is adapted from Press et al.
        !
        l = coefficient_count/2 + 1
        ir = coefficient_count
        do
            if (l > 1) then
                l = l - 1
                cof = coeff(l)
                icl = iclmn(l)
                ind = index(l)
                lab = label(l)
            else
                cof = coeff(ir)
                icl = iclmn(ir)
                ind = index(ir)
                lab = label(ir)
                coeff(ir) = coeff(1)
                iclmn(ir) = iclmn(1)
                index(ir) = index(1)
                label(ir) = label(1)
                ir = ir - 1
                if (ir == 1) then
                    coeff(1) = cof
                    iclmn(1) = icl
                    index(1) = ind
                    label(1) = lab
                    exit
                end if
            end if
            i = l
            j = l + l
            do while (j <= ir)
                if (j < ir) then
                    if (label(j) < label(j+1)) j = j + 1
                end if
                if (lab < label(j)) then
                    coeff(i) = coeff(j)
                    iclmn(i) = iclmn(j)
                    index(i) = index(j)
                    label(i) = label(j)
                    i = j
                    j = j + j
                else
                    j = ir + 1
                end if
            end do
            coeff(i) = cof
            iclmn(i) = icl
            index(i) = ind
            label(i) = lab
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
        last = label(1)
        ibeg = 1
        iend = 1
        distinct_integral_count = 1
        !
        do i = 2, coefficient_count
            if (label(i) == last) then
                iend = iend + 1
            else
                write (mcp_append_unit) last, iend - ibeg + 1
                !GG               WRITE (9999,*) LAST, IEND - IBEG + 1
                write (mcp_append_unit) (iclmn(j),index(j),coeff(j),j=ibeg,iend)
                !cjb-GG        DO J = IBEG,IEND
                !GG               WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
                !cjb-GG        END DO

                !              IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,IEND)
                distinct_integral_count = distinct_integral_count + 1
                last = label(i)
                ibeg = iend + 1
                iend = ibeg
            end if
        end do
        !
        if (ibeg <= coefficient_count) then
            write (mcp_append_unit) last, coefficient_count - ibeg + 1
            !GG            WRITE (9999,*) LAST, NCOEFF - IBEG + 1
            write (mcp_append_unit) (iclmn(j),index(j),coeff(j),j=ibeg,coefficient_count)
            do j=ibeg,coefficient_count
                !GG            WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
            end do

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
            last = label(1)
            ibeg = 1
            iend = 1
            !
            do i = 2, coefficient_count
                if (label(i) == last) iend = iend + 1
                do
                    if (label(i)==last .and. i/=coefficient_count) exit
                    lab = last
                    ncontr = iend - ibeg + 1
                    if (mcp_channel_index == 31) then
                        ia = mod(lab,orbital_label_radix)
                        ib = lab/orbital_label_radix
                        call ensure_debug_output()
                        write (angular_debug_unit, '(" I(",1i2,1a2,",",1i2,1a2,"):")') &
                            orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                            & orbital_symmetry_labels(ib)
                        do j = ibeg, iend
                            call ensure_debug_output()
                            write (angular_debug_unit, '("  T_[",1i2,",",1i4,"] = ",1p,d19.12)') iclmn(j), index(j), coeff(j)
                        end do
                    else
                        k = mcp_channel_index - 32
                        id = mod(lab,orbital_label_radix)
                        lab = lab/orbital_label_radix
                        ib = mod(lab,orbital_label_radix)
                        lab = lab/orbital_label_radix
                        ic = mod(lab,orbital_label_radix)
                        ia = lab/orbital_label_radix
                        call ensure_debug_output()
                        write (angular_debug_unit, '(" R^[(",1i2,")] (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2,1a2,"):")') &
                            k, orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                            & orbital_symmetry_labels(ib), orbital_principal_n(ic), orbital_symmetry_labels&
                            (ic), orbital_principal_n(id), orbital_symmetry_labels(id)
                        do j = ibeg, iend
                            call ensure_debug_output()
                            write (angular_debug_unit, '("  V^[(",1i2,")]_[",1i8,",",1i8,"] = ",1p,d19.12)') &
                                k, iclmn(j), index(j), coeff(j)
                        end do
                    end if
                    last = label(i)
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
        call dalloc (coeff,"COEFF", "SORT")
        call dalloc (iclmn, "ICLMN", "SORT")
        call dalloc (index, "INDEX", "SORT")
        call dalloc (label, "LABEL", "SORT")
        if (mcp_channel_index == 33) call dalloc (nswap, "NSWAP", "SORT")
    end if


    return
end subroutine sort_mcp_coefficients
