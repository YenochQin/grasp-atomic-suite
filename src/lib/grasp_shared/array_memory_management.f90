! Allocation, resizing, and deallocation of pointer targets with caller-specific diagnostics.
module array_memory_management
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64, int64, int8, output_unit
    implicit none (type, external)

    ! Module interface
    private :: real64, int64, int8, output_unit

    ! Procedure interfaces
    ! Supply a pointer with no live target to alloc. Use ralloc when existing data must be kept.
    interface alloc
        module procedure alloc_0i, alloc_1i, alloc_2i, alloc_3i,     &
            alloc_1il,                        &
            alloc_1r, alloc_2r,                         &
            alloc_1rl,                                  &
            alloc_0b, alloc_1b, alloc_2b, alloc_3b,     &
            alloc_1c
    end interface

    ! Procedure interfaces
    ! Resizing preserves the overlapping elements of an associated target.
    interface ralloc
        module procedure ralloc_0i, ralloc_1i, ralloc_2i, ralloc_3i, &
            ralloc_1r, ralloc_2r,                       &
            ralloc_0b, ralloc_1b, ralloc_2b, ralloc_3b, &
            ralloc_1c
    end interface ralloc

    ! Procedure interfaces
    interface dalloc
        module procedure dalloc_1i, dalloc_2i, dalloc_3i,            &
            dalloc_1r, dalloc_2r,                       &
            dalloc_1b, dalloc_2b, dalloc_3b,            &
            dalloc_1c
    end interface

contains
    subroutine alloc_0i(array_pointer, lower_bound, upper_bound, variable_label, caller_label )
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:)
        integer, intent(in) :: lower_bound, upper_bound
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(lower_bound:upper_bound),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,2I15)") "Number of elements requested: ",lower_bound, upper_bound
                stop
            end if
        end if
    end subroutine alloc_0i

    subroutine alloc_1i(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:)
        integer, intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(element_count),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
        end if
    end subroutine alloc_1i

    subroutine alloc_2i(array_pointer, first_dimension_size, second_dimension_size, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(first_dimension_size,second_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15,I15)") "Number of elements requested: ",first_dimension_size,second_dimension_size
                stop
            end if
        end if
    end subroutine alloc_2i

    subroutine alloc_3i(array_pointer, first_dimension_size, second_dimension_size, third_dimension_size, &
        & variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:, :, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size, third_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(first_dimension_size,second_dimension_size,third_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,3I15)") "Number of elements requested: ",first_dimension_size,second_dimension_size, &
                    third_dimension_size
                stop
            end if
        end if
    end subroutine alloc_3i

    subroutine alloc_1il(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:)
        integer(int64), intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(element_count),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
        end if
    end subroutine alloc_1il

    subroutine alloc_1r(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        real(real64), pointer :: array_pointer(:)
        integer, intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(element_count),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
        end if
    end subroutine alloc_1r

    subroutine alloc_2r(array_pointer, first_dimension_size, second_dimension_size, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        real(real64), pointer :: array_pointer(:, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(first_dimension_size,second_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15,I15)") "Number of elements requested: ",first_dimension_size, second_dimension_size
                stop
            end if
        end if
    end subroutine alloc_2r

    subroutine alloc_1rl(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        real(real64), pointer :: array_pointer(:)
        integer(int64), intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(element_count),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
        end if
    end subroutine alloc_1rl

    subroutine alloc_0b(array_pointer, lower_bound, upper_bound, variable_label, caller_label )
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:)
        integer, intent(in) :: lower_bound, upper_bound
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(lower_bound:upper_bound),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,2I15)") "Number of elements requested: ",lower_bound, upper_bound
                stop
            end if
        end if
    end subroutine alloc_0b

    subroutine alloc_1b(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:)
        integer, intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(element_count),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
        end if
    end subroutine alloc_1b

    subroutine alloc_2b(array_pointer, first_dimension_size, second_dimension_size, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(first_dimension_size,second_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15,I15)") "Number of elements requested: ",first_dimension_size,second_dimension_size
                stop
            end if
        end if
    end subroutine alloc_2b

    subroutine alloc_3b(array_pointer, first_dimension_size, second_dimension_size, third_dimension_size, &
        & variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:, :, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size, third_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        nullify(array_pointer)
        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(first_dimension_size,second_dimension_size,third_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,3I15)") "Number of elements requested: ",first_dimension_size,second_dimension_size, &
                    third_dimension_size
                stop
            end if
        end if
    end subroutine alloc_3b

    subroutine alloc_1c(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        character(len=*), pointer :: array_pointer(:)
        integer, intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status

        if (associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " already allocated"
            stop
        else
            allocate(array_pointer(element_count),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
        end if
    end subroutine alloc_1c

    subroutine dalloc_1i(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_1i

    subroutine dalloc_2i(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:, :)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_2i

    subroutine dalloc_3i(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:, :, :)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_3i

    subroutine dalloc_1r(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        real(real64), pointer :: array_pointer(:)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_1r

    subroutine dalloc_2r(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        real(real64), pointer :: array_pointer(:, :)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_2r

    subroutine dalloc_1b(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_1b

    subroutine dalloc_2b(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:, :)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_2b

    subroutine dalloc_3b(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:, :, :)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_3b

    subroutine dalloc_1c(array_pointer, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        character(len=256), pointer :: array_pointer(:)
        character(len=*) :: variable_label, caller_label
        if (associated(array_pointer)) then
            deallocate (array_pointer)
        else
            write(output_unit,"(T10,A,A,A,A,A)") " DALLOC:  Variable ",trim(variable_label), &
                " in subroutine ", trim(caller_label),  " not allocated"
        end if
    end subroutine dalloc_1c


    ! The explicit-bound integer and byte overloads compute the old upper index as SIZE - 1.
    subroutine ralloc_0i(array_pointer, lower_bound, upper_bound, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:), resized_array_pointer(:)
        integer, intent(in) :: lower_bound, upper_bound
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, last_preserved_index
        if (.not.  associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(lower_bound:upper_bound),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to reallocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Elements requested for range: ", &
                    lower_bound, " to ", upper_bound
                stop
            end if
            last_preserved_index = min(size(array_pointer), upper_bound+1) - 1
            resized_array_pointer(lower_bound:last_preserved_index) = array_pointer(lower_bound:last_preserved_index)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_0i


    subroutine ralloc_1i(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:), resized_array_pointer(:)
        integer, intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, preserved_element_count
        if (.not.  associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(element_count),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to reallocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
            preserved_element_count = min(size(array_pointer),element_count)
            resized_array_pointer(1:preserved_element_count) = array_pointer(1:preserved_element_count)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_1i

    subroutine ralloc_2i(array_pointer, first_dimension_size, second_dimension_size, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:, :), resized_array_pointer(:, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, n1_old, n2_old

        if (.not. associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(first_dimension_size,second_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15,I15)") "Number of elements requested: ",first_dimension_size,second_dimension_size
                stop
            end if
            n1_old = min(size(array_pointer, dim=1),first_dimension_size)
            n2_old = min(size(array_pointer, dim=2),second_dimension_size)
            resized_array_pointer(1:n1_old, 1:n2_old) = array_pointer(1:n1_old, 1:n2_old)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_2i

    subroutine ralloc_3i(array_pointer, first_dimension_size, second_dimension_size, third_dimension_size, &
        & variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer, pointer :: array_pointer(:, :, :), resized_array_pointer(:, :, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size, third_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, n1_old, n2_old, n3_old

        if (.not. associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(first_dimension_size,second_dimension_size,third_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15,I15)") "Number of elements requested: ",first_dimension_size,second_dimension_size, &
                    & third_dimension_size
                stop
            end if
            n1_old = min(size(array_pointer, dim=1),first_dimension_size)
            n2_old = min(size(array_pointer, dim=2),second_dimension_size)
            n3_old = min(size(array_pointer, dim=3),third_dimension_size)
            resized_array_pointer(1:n1_old, 1:n2_old, 1:n3_old) = array_pointer(1:n1_old, 1:n2_old, 1:n3_old)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_3i

    subroutine ralloc_1r(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        real(real64), pointer :: array_pointer(:), resized_array_pointer(:)
        integer, intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, preserved_element_count

        if (.not. associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(element_count),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
            preserved_element_count = min(size(array_pointer),element_count)
            resized_array_pointer(1:preserved_element_count) = array_pointer(1:preserved_element_count)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_1r

    subroutine ralloc_2r(array_pointer, first_dimension_size, second_dimension_size, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        real(real64), pointer :: array_pointer(:, :), resized_array_pointer(:, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, n1_old, n2_old

        if (.not. associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(first_dimension_size,second_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15,I15)") "Number of elements requested: ",first_dimension_size, second_dimension_size
                stop
            end if
            n1_old = min(size(array_pointer, dim=1),first_dimension_size)
            n2_old = min(size(array_pointer, dim=2),second_dimension_size)
            resized_array_pointer(1:n1_old, 1:n2_old) = array_pointer(1:n1_old, 1:n2_old)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_2r

    subroutine ralloc_0b(array_pointer, lower_bound, upper_bound, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:), resized_array_pointer(:)
        integer, intent(in) :: lower_bound, upper_bound
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, last_preserved_index
        if (.not.  associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(lower_bound:upper_bound),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to reallocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Elements requested for range: ", &
                    lower_bound, " to ", upper_bound
                stop
            end if
            last_preserved_index = min(size(array_pointer), upper_bound+1) - 1
            resized_array_pointer(lower_bound:last_preserved_index) = array_pointer(lower_bound:last_preserved_index)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_0b

    subroutine ralloc_1b(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:), resized_array_pointer(:)
        integer, intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, preserved_element_count
        if (.not.  associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(element_count),stat=allocation_status)
            if (allocation_status/= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to reallocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
            preserved_element_count = min(size(array_pointer),element_count)
            resized_array_pointer(1:preserved_element_count) = array_pointer(1:preserved_element_count)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_1b

    subroutine ralloc_2b(array_pointer, first_dimension_size, second_dimension_size, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:, :), resized_array_pointer(:, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, n1_old, n2_old

        if (.not. associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(first_dimension_size,second_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15,I15)") "Number of elements requested: ",first_dimension_size,second_dimension_size
                stop
            end if
            n1_old = min(size(array_pointer, dim=1),first_dimension_size)
            n2_old = min(size(array_pointer, dim=2),second_dimension_size)
            resized_array_pointer(1:n1_old, 1:n2_old) = array_pointer(1:n1_old, 1:n2_old)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_2b

    subroutine ralloc_3b(array_pointer, first_dimension_size, second_dimension_size, third_dimension_size, &
        & variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        integer(int8), pointer :: array_pointer(:, :, :), resized_array_pointer(:, :, :)
        integer, intent(in) :: first_dimension_size, second_dimension_size, third_dimension_size
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, n1_old, n2_old, n3_old, i

        if (.not. associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(first_dimension_size,second_dimension_size,third_dimension_size),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15,I15)") "Number of elements requested: ",first_dimension_size,second_dimension_size, &
                    & third_dimension_size
                stop
            end if
            n1_old = min(size(array_pointer, dim=1),first_dimension_size)
            n2_old = min(size(array_pointer, dim=2),second_dimension_size)
            n3_old = min(size(array_pointer, dim=3),third_dimension_size)
            do i = 1,n2_old
                resized_array_pointer(1:n1_old, i, 1:n3_old) = array_pointer(1:n1_old, i, 1:n3_old)
            end do
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_3b

    subroutine ralloc_1c(array_pointer, element_count, variable_label, caller_label)
        implicit none (type, external)

        ! Arguments
        character(len=256), pointer :: array_pointer(:), resized_array_pointer(:)
        integer, intent(in) :: element_count
        character(len=*), intent(in) :: variable_label, caller_label

        ! Local declarations
        integer :: allocation_status, preserved_element_count

        if (.not. associated(array_pointer)) then
            write(output_unit,"(T10,A,A,A,A,A)") "Variable ",trim(variable_label)," in subroutine ", &
                trim(caller_label),  " not already allocated"
            stop
        else
            allocate(resized_array_pointer(element_count),stat=allocation_status)
            if (allocation_status /= 0) then
                write(output_unit,"(T10,A,A,A,A)") "Unable to allocate ",trim(variable_label), &
                    " in subroutine ", trim(caller_label)
                write(output_unit, "(T10,A,I15)") "Number of elements requested: ",element_count
                stop
            end if
            preserved_element_count = min(size(array_pointer),element_count)
            resized_array_pointer(1:preserved_element_count) = array_pointer(1:preserved_element_count)
            deallocate(array_pointer)
            array_pointer => resized_array_pointer
        end if
    end subroutine ralloc_1c

end module array_memory_management
