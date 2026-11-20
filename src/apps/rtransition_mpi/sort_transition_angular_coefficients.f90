!***********************************************************************
!                                                                      *
subroutine sort_transition_angular_coefficients(name, coefficient_scratch_unit, transition_angular_unit, print_details, &
    & transition_type_index, iblki, iblkf)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: iostat_end, output_unit
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   Routine to sort angular coefficients into list based on integral   *
    !   labels rather than CSF.  A tree sort is used. To save space, the   *
    !   CSF pair labels and coefficients  are not read into arrays until   *
    !   the sorting has been done.                                         *
    !                                                                      *
    !   Call(s) to: [OSCL92]: ALCLLA, ALCNMA.                              *
    !                                                                      *
    !   Original author(s) unknown. Modifications for dynamic storage by   *
    !   Farid A. Parpia.                                                   *
    !                                                                      *
    !                                           Last update: 28 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: orbital_index_radix
    use array_memory_management, only: alloc, dalloc, ralloc
    use calculation_input_defaults, only: inpci, ndef, ndump
    use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_symmetry_labels
    use transition_matrix_workspace, only: nkp, kp
    use mpi_state, only: mpi_process_count, mpi_rank

    ! Procedure interfaces
    use resize_transition_label_arrays_interface, only: resize_transition_label_arrays
    use resize_transition_coefficient_arrays_interface, only: resize_transition_coefficient_arrays
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: coefficient_scratch_unit
    integer, intent(inout) :: transition_angular_unit
    integer  :: transition_type_index
    integer, intent(in) :: iblki
    integer, intent(in) :: iblkf
    logical, intent(in) :: print_details
    character, intent(in) :: name(2)*24

    ! Local constants
    integer, parameter :: nca = 65536
    integer, parameter :: orbital_label_radix = orbital_index_radix

    ! Local declarations
    integer, pointer :: ilab(:), iright(:), ileft(:), ibeg(:), iptr(:), &
        ilast(:), iblint(:), iptcsf(:), lblint(:), jlabl(:), isldr(:), isldr1(:)
    real(real64), pointer :: xsldr(:), xl(:)
    integer :: nmcp, nint, nlabel, lldim, nmdim, ir, is, ni, i, m, j, &
        icount, jlab, k, l, imcp, j1, j2, mlr, mup, ints, ia, ib
    integer :: io_status
    character(len=256) :: io_message
    logical :: first
    character(len=2) :: s(-9:9)
    !-----------------------------------------------
    !
    !
    !     POINTER (PIBEG,IBEG(1)),(PILAB,ILAB(1)),(PILAST,ILAST(1)),
    !    :        (PILEFT,ILEFT(1)),(PIPTCS,IPTCSF(1)),(PIRIGH,IRIGHT(1)),
    !    :        (PLBLIN,LBLINT(1))
    !     POINTER (PIPTR,IPTR(1))
    !
    !
    s((-9)) = "-9"
    s((-8)) = "-8"
    s((-7)) = "-7"
    s((-6)) = "-6"
    s((-5)) = "-5"
    s((-4)) = "-4"
    s((-3)) = "-3"
    s((-2)) = "-2"
    s((-1)) = "-1"
    s(0) = "+0"
    s(1) = "+1"
    s(2) = "+2"
    s(3) = "+3"
    s(4) = "+4"
    s(5) = "+5"
    s(6) = "+6"
    s(7) = "+7"
    s(8) = "+8"
    s(9) = "+9"
    !
    !   Position file at beginning of list of integrals
    !
    rewind (coefficient_scratch_unit)
    !
    !   Initialize
    !
    first = .true.
    nmcp = 0
    nint = 0
    !
    !   Initial allocation of storage to local arrays
    !
    nlabel = 1
    call alloc (jlabl, nlabel, "JLABL", "TRSORT")
    call alloc (xl, nlabel, "XL", "TRSORT")
    !
    call resize_transition_label_arrays (ibeg, ilab, ilast, ileft, iptcsf, iright, lblint, &
        lldim, 1)
    call resize_transition_coefficient_arrays (iptr, isldr, isldr1, xsldr, nmdim, 1)
    !
    !   Now the rest of the elements
    !
    m = 0
    do
        read (coefficient_scratch_unit, iostat=io_status, iomsg=io_message) ir, is, ni
        if (io_status == iostat_end) exit
        if (io_status /= 0) error stop trim(io_message)
        if (ni > nlabel) then
            call ralloc (jlabl, ni, "JLABL", "TRSORT")
            call ralloc (xl, ni, "XL", "TRSORT")
            nlabel = ni
        end if
        read (coefficient_scratch_unit, iostat=io_status) (jlabl(i),xl(i),i=1,ni)
        if (io_status /= 0) then
            write (output_unit, *) "TRSORT: Error reading CSF-based file."
            stop
        end if
        if (ir==0 .or. is==0 .or. ni==0) cycle

        icount = 0
        if (first) then
            ! Find the first nonzero element and initialize the tree.
            do icount = 1, ni
                if (jlabl(icount) /= 0) exit
            end do
            if (icount > ni) cycle
            ilab(1) = jlabl(icount)
            iright(1) = 0
            ileft(1) = 0
            ibeg(1) = 1
            iptr(1) = 0
            ilast(1) = 0
            m = 1
            j = 1
            first = .false.
        end if

        ! Add each remaining nonzero element to the integral tree.
        do icount = icount + 1, ni
            jlab = jlabl(icount)
            if (jlab == 0) cycle

            m = m + 1
            if (m > nmdim) call resize_transition_coefficient_arrays (iptr, isldr, isldr1, xsldr, nmdim, 2)
            i = 1
            do while (jlab /= ilab(i))
                if (jlab > ilab(i)) then
                    k = ileft(i)
                else
                    k = iright(i)
                end if
                if (k == 0) exit
                i = k
            end do

            if (jlab == ilab(i)) then
                ! Append this CSF pair to the existing integral.
                l = ibeg(i)
                do while (iptr(l) /= 0)
                    l = iptr(l)
                end do
                iptr(l) = m
            else
                j = j + 1
                if (j > lldim) call resize_transition_label_arrays (ibeg, ilab, ilast, ileft, iptcsf, &
                    iright, lblint, lldim, 2)
                if (jlab > ilab(i)) then
                    ileft(i) = j
                else
                    iright(i) = j
                end if
                ilast(j) = i
                iright(j) = 0
                ileft(j) = 0
                ibeg(j) = m
                ilab(j) = jlab
            end if
            iptr(m) = 0
        end do
    end do

    if (.not. first) then
        ! Unpack the tree in ascending order of integral label.
        nmcp = m
        nint = j
        l = 0
        m = 0
        i = 1
        do while (iright(i) /= 0)
            i = iright(i)
        end do
        do while (i /= 0)
            if (ilab(i) == 0) then
                i = ilast(i)
                cycle
            end if
            l = l + 1
            lblint(l) = ilab(i)
            k = ibeg(i)
            do while (k /= 0)
                m = m + 1
                isldr(m) = k
                k = iptr(k)
            end do
            iptcsf(l) = m
            ilab(i) = 0

            if (ileft(i) == 0) then
                i = ilast(i)
            else
                i = ileft(i)
                do while (iright(i) /= 0)
                    i = iright(i)
                end do
            end if
        end do

        ! Invert the CSF/coefficient pointers to obtain their sorted positions.
        do i = 1, nmcp
            k = isldr(i)
            iptr(k) = i
        end do

        ! Read the CSF pairs and coefficients into the sorted list.
        imcp = 0
        rewind (coefficient_scratch_unit)
        do
            read (coefficient_scratch_unit, iostat=io_status, iomsg=io_message) ir, is, ni
            if (io_status == iostat_end) exit
            if (io_status /= 0) error stop trim(io_message)
            if (ir==0 .or. is==0 .or. ni==0) cycle
            read (coefficient_scratch_unit) (jlabl(i),xl(i),i=1,ni)
            do i = 1, ni
                if (jlabl(i) == 0) cycle
                imcp = imcp + 1
                k = iptr(imcp)
                isldr(k) = ir
                isldr1(k) = is
                !              ISLDR(K) = IR*NCA+IS
                xsldr(k) = xl(i)
            end do
        end do
    end if
    rewind (coefficient_scratch_unit)
    !
    !  If first set of data open the file and print
    !  some data to later be able to identify the file
    !
    if (iblki==1 .and. iblkf==1) then
        j1 = index(name(1)," ")
        j2 = index(name(2)," ")
        open(newunit=transition_angular_unit,file=name(1)(1:j1-1)//"."//name(2)(1:j2-1)//"."//s(&
            kp(transition_type_index))//"T", status="UNKNOWN", form="UNFORMATTED", position=&
            "asis")
    end if
    write (transition_angular_unit) iblki,iblkf,orbital_count,nkp,mpi_process_count,mpi_rank
    write (transition_angular_unit) nint
    if (nmcp /= 0) then
        if (print_details) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,/,"  k"/," d  (rs)  Coefficients:"/,"  ab"/,/,/,&
                &"    a     b          r         s  Coefficient"/)')
        end if
        mlr = 1
        do i = 1, nint
            mup = iptcsf(i)
            ints = mup - mlr + 1
            write (transition_angular_unit) lblint(i), ints
            if (print_details) then
                ia = lblint(i)/orbital_label_radix
                ib = mod(lblint(i),orbital_label_radix)
                call ensure_debug_output()
                write (angular_debug_unit, '(2(2x,i2,a2))') &
                    orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                    & orbital_symmetry_labels(ib)
            end if
            write (transition_angular_unit) (isldr(m),isldr1(m),xsldr(m),m=mlr,mup)
            if (print_details) then
                do m = mlr, mup
                    !                  IS = MOD (ISLDR(M),NCA)
                    is = isldr1(m)
                    !                 IR = ISLDR(M)/NCA
                    ir = isldr(m)
                    call ensure_debug_output()
                    write (angular_debug_unit, '(14x,1i6,2x,1i6,2x,1p,1d22.15)') ir, is, xsldr(m)
                end do
            end if
            mlr = mup + 1
        end do
        !
    end if
    !
    !   Deallocate storage for local arrays
    !
    call dalloc (jlabl, "JLABJ", "TRSORT")
    call dalloc (xl, "XL", "TRSORT")
    call resize_transition_label_arrays (ibeg, ilab, ilast, ileft, iptcsf, iright, lblint, &
        lldim, 3)
    call resize_transition_coefficient_arrays (iptr, isldr, isldr1, xsldr, nmdim, 3)
    !
    return
    !
end subroutine sort_transition_angular_coefficients
