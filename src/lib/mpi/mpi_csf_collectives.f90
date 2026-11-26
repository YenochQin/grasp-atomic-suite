! Synchronize CSF descriptor errors and broadcast arrays in bounded batches.
module mpi_csf_collectives
    ! Module dependencies
    use iso_fortran_env, only: int8
    use mpi_f08
    use mpi_element_counts, only: max_mpi_count
    use grasp_parameters, only: max_orbitals
    use csf_coupling_data, only: csf_angular_code_count
    implicit none (type, external)

    ! Module interface
    private
    public :: require_collective_success, allocate_descriptor_arrays, broadcast_descriptor_arrays
contains
    subroutine require_collective_success(status)
        ! Arguments
        integer, intent(in) :: status

        ! Local declarations
        integer :: rank_failed, any_rank_failed, mpi_return_code
        ! Reduce failure flags so every rank takes the same error path.
        rank_failed=merge(1,0,status /= 0)
        call MPI_Allreduce(rank_failed,any_rank_failed,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_return_code)
        if (any_rank_failed /= 0 .or. mpi_return_code /= MPI_SUCCESS) then
            call MPI_Abort(MPI_COMM_WORLD,1,mpi_return_code)
            error stop "Descriptor collective failed"
        end if
    end subroutine require_collective_success

    subroutine allocate_descriptor_arrays(csf_count,subshell_occupations,subshell_quantum_codes,angular_codes)
        ! Arguments
        integer, intent(in) :: csf_count
        integer(int8), allocatable, intent(out) :: subshell_occupations(:,:),subshell_quantum_codes(:,:,:),angular_codes(:,:)

        ! Local declarations
        integer :: allocation_status
        call require_collective_success(merge(0,1,csf_count > 0))
        allocate(subshell_occupations(max_orbitals,csf_count),stat=allocation_status)
        if (allocation_status == 0) allocate(subshell_quantum_codes(max_orbitals,3,csf_count),stat=allocation_status)
        if (allocation_status == 0) allocate(angular_codes(csf_angular_code_count,csf_count),stat=allocation_status)
        call require_collective_success(allocation_status)
    end subroutine allocate_descriptor_arrays

    subroutine broadcast_descriptor_arrays(subshell_occupations,subshell_quantum_codes,angular_codes,csf_batch_limit)
        ! Arguments
        ! Fortran packs strided actual arguments on entry. MPI sees contiguous batches.
        integer(int8), contiguous, intent(inout) :: subshell_occupations(:,:)
        integer(int8), contiguous, intent(inout), optional :: subshell_quantum_codes(:,:,:),angular_codes(:,:)
        integer, intent(in), optional :: csf_batch_limit

        ! Local declarations
        integer :: local_array_layout(4), root_array_layout(4), layout_status
        integer :: first_csf_index, last_csf_index, csfs_per_batch, mpi_return_code
        integer(mpi_count_kind) :: max_elements_per_csf
        ! Ranks must agree on which optional arrays are present before entering the payload collectives.
        local_array_layout = [size(subshell_occupations,1), size(subshell_occupations,2), &
            merge(1,0,present(subshell_quantum_codes)), merge(1,0,present(angular_codes))]
        root_array_layout=local_array_layout
        call MPI_Bcast(root_array_layout,4,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        layout_status=merge(0,1,all(root_array_layout == local_array_layout))
        if (present(subshell_quantum_codes)) then
            if (any(shape(subshell_quantum_codes) /= [local_array_layout(1),3,local_array_layout(2)])) layout_status=1
        end if
        if (present(angular_codes)) then
            if (any(shape(angular_codes) /= [csf_angular_code_count,local_array_layout(2)])) layout_status=1
        end if
        if (local_array_layout(1) <= 0) layout_status=1
        call require_collective_success(layout_status)
        ! Limit the MPI element count using the widest array per CSF.
        max_elements_per_csf=int(local_array_layout(1),mpi_count_kind)
        if (present(subshell_quantum_codes)) max_elements_per_csf=3_mpi_count_kind*max_elements_per_csf
        if (present(angular_codes)) max_elements_per_csf=max(max_elements_per_csf,int(csf_angular_code_count,mpi_count_kind))
        csfs_per_batch=int(min(8192_mpi_count_kind,max_mpi_count/max_elements_per_csf))
        if (present(csf_batch_limit)) csfs_per_batch=min(csfs_per_batch,csf_batch_limit)
        call require_collective_success(merge(0,1,csfs_per_batch > 0))
        ! Require every rank to use the same sequence of batch sizes.
        root_array_layout(1)=csfs_per_batch
        call MPI_Bcast(root_array_layout,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_return_code)
        call require_collective_success(merge(0,1,csfs_per_batch == root_array_layout(1)))
        ! Slice batches along the last dimension so each transfer contains complete CSF records.
        first_csf_index=1
        do while (first_csf_index <= local_array_layout(2))
            last_csf_index=int(min(int(local_array_layout(2),mpi_count_kind),int(first_csf_index,mpi_count_kind)+csfs_per_batch-1))
            call MPI_Bcast(subshell_occupations(:,first_csf_index:last_csf_index), &
                local_array_layout(1)*(last_csf_index-first_csf_index+1), MPI_INTEGER1, 0, MPI_COMM_WORLD, mpi_return_code)
            if (mpi_return_code /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,1,mpi_return_code)
            if (present(subshell_quantum_codes)) then
                call MPI_Bcast(subshell_quantum_codes(:,:,first_csf_index:last_csf_index), &
                    3*local_array_layout(1)*(last_csf_index-first_csf_index+1), &
                    MPI_INTEGER1,0,MPI_COMM_WORLD,mpi_return_code)
                if (mpi_return_code /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,1,mpi_return_code)
            end if
            if (present(angular_codes)) then
                call MPI_Bcast(angular_codes(:,first_csf_index:last_csf_index), &
                    csf_angular_code_count*(last_csf_index-first_csf_index+1), &
                    MPI_INTEGER1, 0, MPI_COMM_WORLD, mpi_return_code)
                if (mpi_return_code /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD,1,mpi_return_code)
            end if
            if (last_csf_index == local_array_layout(2)) exit
            first_csf_index=last_csf_index+1
        end do
    end subroutine broadcast_descriptor_arrays
end module mpi_csf_collectives
