! Merge two CSF bases into a common orbital order, with the left basis preceding the right.
module csf_basis_merging
    ! Module dependencies
    use iso_fortran_env, only: int8, int64
    use grasp_parameters, only: max_orbitals
    use atomic_calculation_state, only: electron_count, left_electron_count
    use orbital_csf_basis, only: csf_count, left_csf_count, left_orbital_count, left_principal_n, left_kappa, &
        left_orbital_labels, left_occupations, orbital_count, orbital_principal_n, orbital_kappa, orbital_l, orbital_two_j, &
        orbital_symmetry_labels, csf_subshell_occupations
    use csf_coupling_data, only: left_subshell_quantum_codes, left_angular_codes, csf_subshell_quantum_codes, &
        csf_angular_codes, csf_angular_code_count, csf_total_j_slot
    implicit none (type, external)

    ! Module interface
    private
    public :: merge_csf_bases
contains
    subroutine merge_csf_bases(left_core_subshell_count,right_core_subshell_count,status,message)
        ! Arguments
        integer, intent(in) :: left_core_subshell_count,right_core_subshell_count
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        ! Local declarations
        integer :: left_orbital_map(max_orbitals), right_orbital_map(max_orbitals)
        integer :: merged_principal_n(max_orbitals), merged_kappa(max_orbitals)
        character(len=2) :: merged_labels(max_orbitals)
        integer :: left_orbital_index, right_orbital_index, merge_choice, merged_orbital_count, merged_orbital_index, &
            & input_orbital_index, merged_csf_count, allocation_status
        integer(int64) :: wide_merged_csf_count
        logical :: left_orbital_in_right_tail, right_orbital_in_left_tail
        integer(int8), allocatable :: occupations(:,:), subshell_quantum_codes(:,:,:), angular_codes(:,:)
        status=1
        message="left and right electron counts or cores differ"
        if (electron_count /= left_electron_count .or. left_core_subshell_count /= right_core_subshell_count) return
        if (left_core_subshell_count < 0 .or. left_core_subshell_count > min(left_orbital_count,orbital_count)) return
        if (any(left_principal_n(:left_core_subshell_count) /= orbital_principal_n(:right_core_subshell_count)) .or. &
            any(left_kappa(:left_core_subshell_count) /= orbital_kappa(:right_core_subshell_count))) return
        message="input basis arrays are not loaded"
        if (.not.allocated(left_occupations) .or. .not.allocated(left_subshell_quantum_codes) .or. &
            & .not.allocated(left_angular_codes)) return
        if (.not.allocated(csf_subshell_occupations) .or. .not.allocated(csf_subshell_quantum_codes) .or. &
            .not.allocated(csf_angular_codes)) return

        left_orbital_index=1
        right_orbital_index=1
        merged_orbital_count=0
        ! Merge the ordered orbital lists. Prefer the left orbital if neither order constrains the choice.
        do while (left_orbital_index <= left_orbital_count .or. right_orbital_index <= orbital_count)
            if (left_orbital_index > left_orbital_count) then
                merge_choice=2
            else if (right_orbital_index > orbital_count) then
                merge_choice=1
            else if (left_principal_n(left_orbital_index) == orbital_principal_n(right_orbital_index) .and. &
                & left_kappa(left_orbital_index) == orbital_kappa(right_orbital_index)) then
                merge_choice=3
            else
                left_orbital_in_right_tail=any(left_principal_n(left_orbital_index) == &
                    & orbital_principal_n(right_orbital_index:orbital_count) .and. &
                    left_kappa(left_orbital_index) == orbital_kappa(right_orbital_index:orbital_count))
                right_orbital_in_left_tail=any(orbital_principal_n(right_orbital_index) == &
                    & left_principal_n(left_orbital_index:left_orbital_count) .and. &
                    orbital_kappa(right_orbital_index) == left_kappa(left_orbital_index:left_orbital_count))
                if (.not.left_orbital_in_right_tail) then
                    merge_choice=1
                else if (.not.right_orbital_in_left_tail) then
                    merge_choice=2
                else
                    message="inconsistent orbital order between the two bases"
                    return
                end if
            end if
            message="merged orbital count exceeds configured capacity"
            if (merged_orbital_count == max_orbitals) return
            merged_orbital_count=merged_orbital_count+1
            if (merge_choice == 2) then
                merged_principal_n(merged_orbital_count)=orbital_principal_n(right_orbital_index)
                merged_kappa(merged_orbital_count)=orbital_kappa(right_orbital_index)
                merged_labels(merged_orbital_count)=orbital_symmetry_labels(right_orbital_index)
            else
                merged_principal_n(merged_orbital_count)=left_principal_n(left_orbital_index)
                merged_kappa(merged_orbital_count)=left_kappa(left_orbital_index)
                merged_labels(merged_orbital_count)=left_orbital_labels(left_orbital_index)
            end if
            left_orbital_map(merged_orbital_count)=0
            right_orbital_map(merged_orbital_count)=0
            if (merge_choice /= 2) then
                left_orbital_map(merged_orbital_count)=left_orbital_index
                left_orbital_index=left_orbital_index+1
            end if
            if (merge_choice /= 1) then
                right_orbital_map(merged_orbital_count)=right_orbital_index
                right_orbital_index=right_orbital_index+1
            end if
        end do
        wide_merged_csf_count=int(left_csf_count,int64)+int(csf_count,int64)
        message="merged CSF count exceeds suite integer capacity"
        if (wide_merged_csf_count > int(huge(merged_csf_count),int64) .or. wide_merged_csf_count < 1) return
        merged_csf_count=int(wide_merged_csf_count)
        message="cannot allocate merged CSF arrays"
        allocate(occupations(max_orbitals,merged_csf_count), &
            subshell_quantum_codes(max_orbitals,3,merged_csf_count), &
            angular_codes(csf_angular_code_count,merged_csf_count), stat=allocation_status)
        if (allocation_status /= 0) return
        ! When a side has no orbital at this position it inherits its own preceding 2K+1.
        occupations=0_int8
        subshell_quantum_codes=0_int8
        subshell_quantum_codes(:merged_orbital_count,3,:)=1_int8
        angular_codes=0_int8
        angular_codes(csf_total_j_slot:,:left_csf_count)=int(left_angular_codes(csf_total_j_slot:,:left_csf_count),int8)
        angular_codes(csf_total_j_slot:,left_csf_count+1:)=csf_angular_codes(csf_total_j_slot:,:)
        do merged_orbital_index=1,merged_orbital_count
            if (merged_orbital_index == 1) then
                angular_codes(merged_orbital_index,:)=1_int8
            else
                angular_codes(merged_orbital_index,:)=angular_codes(merged_orbital_index-1,:)
            end if
            input_orbital_index=left_orbital_map(merged_orbital_index)
            if (input_orbital_index /= 0) then
                occupations(merged_orbital_index,:left_csf_count)=int(left_occupations(input_orbital_index,:left_csf_count),int8)
                subshell_quantum_codes(merged_orbital_index,:, &
                    :left_csf_count)=int(left_subshell_quantum_codes(input_orbital_index,:,:left_csf_count),int8)
                angular_codes(merged_orbital_index,:left_csf_count) = &
                    int(left_angular_codes(input_orbital_index,:left_csf_count),int8)
            end if
            input_orbital_index=right_orbital_map(merged_orbital_index)
            if (input_orbital_index /= 0) then
                occupations(merged_orbital_index,left_csf_count+1:)=csf_subshell_occupations(input_orbital_index,:)
                subshell_quantum_codes(merged_orbital_index,:,left_csf_count+1:)=csf_subshell_quantum_codes(input_orbital_index,:,:)
                angular_codes(merged_orbital_index,left_csf_count+1:)=csf_angular_codes(input_orbital_index,:)
            end if
        end do
        call move_alloc(occupations, csf_subshell_occupations)
        call move_alloc(subshell_quantum_codes, csf_subshell_quantum_codes)
        call move_alloc(angular_codes, csf_angular_codes)
        orbital_count=merged_orbital_count
        csf_count=merged_csf_count
        orbital_principal_n(:merged_orbital_count)=merged_principal_n(:merged_orbital_count)
        orbital_kappa(:merged_orbital_count)=merged_kappa(:merged_orbital_count)
        orbital_symmetry_labels(:merged_orbital_count)=merged_labels(:merged_orbital_count)
        orbital_two_j(:merged_orbital_count)=2*abs(merged_kappa(:merged_orbital_count))-1
        where (merged_kappa(:merged_orbital_count) > 0)
            orbital_l(:merged_orbital_count)=merged_kappa(:merged_orbital_count)
        elsewhere
            orbital_l(:merged_orbital_count)=-merged_kappa(:merged_orbital_count)-1
        end where
        status=0
        message=""
    end subroutine merge_csf_bases
end module csf_basis_merging
