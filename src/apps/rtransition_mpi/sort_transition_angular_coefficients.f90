module sort_transition_angular_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: sort_transition_angular_coefficients
contains

    ! Routine to sort angular coefficients into list based on integral
    ! labels rather than CSF.  A tree sort is used. To save space, the
    ! CSF pair labels and coefficients  are not read into arrays until
    ! the sorting has been done.
    !
    ! Original author(s) unknown. Modifications for dynamic storage by
    ! Farid A. Parpia.
    !
    ! Last update: 28 Dec 1992
    subroutine sort_transition_angular_coefficients(name, coefficient_scratch_unit, transition_angular_unit, print_details, &
        & transition_type_index, initial_block_index, final_block_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: iostat_end, output_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use array_memory_management, only: alloc, dalloc, ralloc
        use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_symmetry_labels
        use transition_matrix_workspace, only: requested_multipole_count, requested_multipole_codes
        use mpi_state, only: mpi_process_count, mpi_rank

        ! Procedure interfaces
        use resize_transition_label_arrays_interface, only: resize_transition_label_arrays
        use resize_transition_coefficient_arrays_interface, only: resize_transition_coefficient_arrays
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: coefficient_scratch_unit
        integer, intent(inout) :: transition_angular_unit
        integer  :: transition_type_index
        integer, intent(in) :: initial_block_index
        integer, intent(in) :: final_block_index
        logical, intent(in) :: print_details
        character, intent(in) :: name(2)*24

        ! Local constants
        integer, parameter :: nca = 65536

        ! Local declarations
        integer, pointer :: node_integral_labels(:), right_child_nodes(:), left_child_nodes(:), &
            & first_coefficient_indices(:), coefficient_sort_links(:), parent_nodes(:), coefficient_end_offsets(:), &
            & sorted_integral_labels(:), &
            record_integral_labels(:), initial_csf_index_workspace(:), final_csf_indices(:)
        real(real64), pointer :: sorted_angular_coefficients(:), record_coefficients(:)
        integer :: total_coefficient_count, distinct_integral_count, record_capacity, label_capacity, coefficient_capacity, &
            & initial_csf_index, final_csf_index, record_coefficient_count, i, m, j, &
            record_coefficient_index, packed_integral_label, k, l, original_coefficient_index, initial_name_end, &
                & final_name_end, first_coefficient_index, last_coefficient_index, integral_coefficient_count, ia, ib
        integer :: io_status
        character(len=256) :: io_message
        logical :: first
        character(len=2) :: s(-9:9)
        ! :        (PILEFT,ILEFT(1)),(PIPTCS,IPTCSF(1)),(PIRIGH,IRIGHT(1)),
        ! :        (PLBLIN,LBLINT(1))
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
        ! Position file at beginning of list of integrals
        rewind (coefficient_scratch_unit)
        ! Initialize
        first = .true.
        total_coefficient_count = 0
        distinct_integral_count = 0
        ! Initial allocation of storage to local arrays
        record_capacity = 1
        call alloc (record_integral_labels, record_capacity, "record_integral_labels", "sort_transition_angular_coefficients")
        call alloc (record_coefficients, record_capacity, "record_coefficients", "sort_transition_angular_coefficients")
        call resize_transition_label_arrays (first_coefficient_indices, node_integral_labels, parent_nodes, &
            & left_child_nodes, coefficient_end_offsets, right_child_nodes, sorted_integral_labels, &
            label_capacity, 1)
        call resize_transition_coefficient_arrays (coefficient_sort_links, initial_csf_index_workspace, final_csf_indices, &
            & sorted_angular_coefficients, coefficient_capacity, 1)
        ! Now the rest of the elements
        m = 0
        do
            read (coefficient_scratch_unit, iostat=io_status, iomsg=io_message) initial_csf_index, final_csf_index, &
                & record_coefficient_count
            if (io_status == iostat_end) exit
            if (io_status /= 0) error stop trim(io_message)
            if (record_coefficient_count > record_capacity) then
                call ralloc (record_integral_labels, record_coefficient_count, "record_integral_labels", &
                    "sort_transition_angular_coefficients")
                call ralloc (record_coefficients, record_coefficient_count, "record_coefficients", &
                    "sort_transition_angular_coefficients")
                record_capacity = record_coefficient_count
            end if
            read (coefficient_scratch_unit, iostat=io_status) (record_integral_labels(i),record_coefficients(i),i=1, &
                & record_coefficient_count)
            if (io_status /= 0) then
                write (output_unit, *) "sort_transition_angular_coefficients: Error reading CSF-based file."
                stop
            end if
            if (initial_csf_index==0 .or. final_csf_index==0 .or. record_coefficient_count==0) cycle

            record_coefficient_index = 0
            if (first) then
                ! Find the first nonzero element and initialize the tree.
                do record_coefficient_index = 1, record_coefficient_count
                    if (record_integral_labels(record_coefficient_index) /= 0) exit
                end do
                if (record_coefficient_index > record_coefficient_count) cycle
                node_integral_labels(1) = record_integral_labels(record_coefficient_index)
                right_child_nodes(1) = 0
                left_child_nodes(1) = 0
                first_coefficient_indices(1) = 1
                coefficient_sort_links(1) = 0
                parent_nodes(1) = 0
                m = 1
                j = 1
                first = .false.
            end if

            ! Add each remaining nonzero element to the integral tree.
            do record_coefficient_index = record_coefficient_index + 1, record_coefficient_count
                packed_integral_label = record_integral_labels(record_coefficient_index)
                if (packed_integral_label == 0) cycle

                m = m + 1
                if (m > coefficient_capacity) call resize_transition_coefficient_arrays (coefficient_sort_links, &
                    & initial_csf_index_workspace, final_csf_indices, sorted_angular_coefficients, coefficient_capacity, 2)
                i = 1
                do while (packed_integral_label /= node_integral_labels(i))
                    if (packed_integral_label > node_integral_labels(i)) then
                        k = left_child_nodes(i)
                    else
                        k = right_child_nodes(i)
                    end if
                    if (k == 0) exit
                    i = k
                end do

                if (packed_integral_label == node_integral_labels(i)) then
                    ! Append this CSF pair to the existing integral.
                    l = first_coefficient_indices(i)
                    do while (coefficient_sort_links(l) /= 0)
                        l = coefficient_sort_links(l)
                    end do
                    coefficient_sort_links(l) = m
                else
                    j = j + 1
                    if (j > label_capacity) call resize_transition_label_arrays (first_coefficient_indices, &
                        & node_integral_labels, parent_nodes, left_child_nodes, coefficient_end_offsets, &
                        right_child_nodes, sorted_integral_labels, label_capacity, 2)
                    if (packed_integral_label > node_integral_labels(i)) then
                        left_child_nodes(i) = j
                    else
                        right_child_nodes(i) = j
                    end if
                    parent_nodes(j) = i
                    right_child_nodes(j) = 0
                    left_child_nodes(j) = 0
                    first_coefficient_indices(j) = m
                    node_integral_labels(j) = packed_integral_label
                end if
                coefficient_sort_links(m) = 0
            end do
        end do

        if (.not. first) then
            ! Unpack the tree in ascending order of integral label.
            total_coefficient_count = m
            distinct_integral_count = j
            l = 0
            m = 0
            i = 1
            do while (right_child_nodes(i) /= 0)
                i = right_child_nodes(i)
            end do
            do while (i /= 0)
                if (node_integral_labels(i) == 0) then
                    i = parent_nodes(i)
                    cycle
                end if
                l = l + 1
                sorted_integral_labels(l) = node_integral_labels(i)
                k = first_coefficient_indices(i)
                do while (k /= 0)
                    m = m + 1
                    initial_csf_index_workspace(m) = k
                    k = coefficient_sort_links(k)
                end do
                coefficient_end_offsets(l) = m
                node_integral_labels(i) = 0

                if (left_child_nodes(i) == 0) then
                    i = parent_nodes(i)
                else
                    i = left_child_nodes(i)
                    do while (right_child_nodes(i) /= 0)
                        i = right_child_nodes(i)
                    end do
                end if
            end do

            ! Invert the CSF/coefficient pointers to obtain their sorted positions.
            do i = 1, total_coefficient_count
                k = initial_csf_index_workspace(i)
                coefficient_sort_links(k) = i
            end do

            ! Read the CSF pairs and coefficients into the sorted list.
            original_coefficient_index = 0
            rewind (coefficient_scratch_unit)
            do
                read (coefficient_scratch_unit, iostat=io_status, iomsg=io_message) initial_csf_index, final_csf_index, &
                    & record_coefficient_count
                if (io_status == iostat_end) exit
                if (io_status /= 0) error stop trim(io_message)
                if (initial_csf_index==0 .or. final_csf_index==0 .or. record_coefficient_count==0) cycle
                read (coefficient_scratch_unit) (record_integral_labels(i),record_coefficients(i),i=1,record_coefficient_count)
                do i = 1, record_coefficient_count
                    if (record_integral_labels(i) == 0) cycle
                    original_coefficient_index = original_coefficient_index + 1
                    k = coefficient_sort_links(original_coefficient_index)
                    initial_csf_index_workspace(k) = initial_csf_index
                    final_csf_indices(k) = final_csf_index
                    ! ISLDR(K) = IR*NCA+IS
                    sorted_angular_coefficients(k) = record_coefficients(i)
                end do
            end do
        end if
        rewind (coefficient_scratch_unit)
        ! If first set of data open the file and print
        ! some data to later be able to identify the file
        if (initial_block_index==1 .and. final_block_index==1) then
            initial_name_end = index(name(1)," ")
            final_name_end = index(name(2)," ")
            open(newunit=transition_angular_unit,file=name(1)(1:initial_name_end-1)//"."//name(2)(1:final_name_end-1)//"."//s(&
                requested_multipole_codes(transition_type_index))//"T", status="UNKNOWN", form="UNFORMATTED", position=&
                "asis")
        end if
        write (transition_angular_unit) initial_block_index,final_block_index,orbital_count,requested_multipole_count, &
            mpi_process_count,mpi_rank
        write (transition_angular_unit) distinct_integral_count
        if (total_coefficient_count /= 0) then
            if (print_details) then
                call ensure_debug_output()
                write (angular_debug_unit, '(/,/,"  k"/," d  (rs)  Coefficients:"/,"  ab"/,/,/,&
                    &"    a     b          r         s  Coefficient"/)')
            end if
            first_coefficient_index = 1
            do i = 1, distinct_integral_count
                last_coefficient_index = coefficient_end_offsets(i)
                integral_coefficient_count = last_coefficient_index - first_coefficient_index + 1
                write (transition_angular_unit) sorted_integral_labels(i), integral_coefficient_count
                if (print_details) then
                    ia = sorted_integral_labels(i)/orbital_index_radix
                    ib = mod(sorted_integral_labels(i),orbital_index_radix)
                    call ensure_debug_output()
                    write (angular_debug_unit, '(2(2x,i2,a2))') &
                        orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                        & orbital_symmetry_labels(ib)
                end if
                write (transition_angular_unit) &
                    & (initial_csf_index_workspace(m),final_csf_indices(m),sorted_angular_coefficients(m),m= &
                        & first_coefficient_index,last_coefficient_index)
                if (print_details) then
                    do m = first_coefficient_index, last_coefficient_index
                        ! IS = MOD (ISLDR(M),NCA)
                        final_csf_index = final_csf_indices(m)
                        ! IR = ISLDR(M)/NCA
                        initial_csf_index = initial_csf_index_workspace(m)
                        call ensure_debug_output()
                        write (angular_debug_unit, '(14x,1i6,2x,1i6,2x,1p,1d22.15)') initial_csf_index, final_csf_index, &
                            & sorted_angular_coefficients(m)
                    end do
                end if
                first_coefficient_index = last_coefficient_index + 1
            end do
        end if
        ! Deallocate storage for local arrays
        call dalloc (record_integral_labels, "record_integral_labels", "sort_transition_angular_coefficients")
        call dalloc (record_coefficients, "record_coefficients", "sort_transition_angular_coefficients")
        call resize_transition_label_arrays (first_coefficient_indices, node_integral_labels, parent_nodes, &
            & left_child_nodes, coefficient_end_offsets, right_child_nodes, sorted_integral_labels, &
            label_capacity, 3)
        call resize_transition_coefficient_arrays (coefficient_sort_links, initial_csf_index_workspace, final_csf_indices, &
            & sorted_angular_coefficients, coefficient_capacity, 3)
        return
    end subroutine sort_transition_angular_coefficients
end module sort_transition_angular_coefficients_interface
