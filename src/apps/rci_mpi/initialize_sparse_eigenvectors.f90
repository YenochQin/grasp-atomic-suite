module initialize_sparse_eigenvectors_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_sparse_eigenvectors
contains

    !***********************************************************************
    !
    subroutine initialize_sparse_eigenvectors (nmax, csf_count, mpi_rank, mpi_process_count,                   &
        initial_vector_count, basis, restart_matrix_unit, hamiltonian_energy_shift)
        ! Module dependencies
        !
        !    Routine for providing initial estimates from upper left corner
        !       of the matrix. It is exact (not estimates) if ncf <= nmax which
        !    is currently set to be 400 in the calling routine.
        !
        !    Matrix is sparse and on the disk
        !
        !   Block version by Xinghong He            Last revision: 14 Dec 1998
        !
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc

        ! Procedure interfaces
        use blas_lapack_interfaces, only: dcopy, dspevx

        ! Module dependencies
        use real_vector_initialization, only: initialize_real_vector
        use mpi_utilities, only: allreduce_mpi_real_sum

        ! Procedure interfaces
        use blas_lapack_interfaces, only: slamch
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nmax, csf_count, mpi_rank, mpi_process_count, initial_vector_count, restart_matrix_unit
        real(real64), intent(in) :: hamiltonian_energy_shift
        real(real64), intent(inout):: basis(*)

        ! Local declarations
        real(real64), pointer :: ap(:), eigval(:), vec(:), real_workspace(:), hmx(:)
        integer, pointer :: integer_workspace(:), ifail(:), irow(:)
        integer :: in, ns, file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count, info, &
            ierr, j, joff, nelc, ir, nfound
        real(real64) :: elsto, abstol
        !-----------------------------------------------------------------------
        ns = min (nmax, csf_count)

        call alloc (ap, (ns*(ns+1))/2, "AP", "initialize_sparse_eigenvectors")
        !GG      CALL alloc (ap_buffer, (NS*(NS+1))/2, 'AP', 'INIESTSD')
        call initialize_real_vector ((ns*(ns+1))/2, 0.0_real64, ap, 1)

        !**** separate upper left block of size NS*NS

        call alloc (hmx, csf_count, "HMX", "initialize_sparse_eigenvectors")
        call alloc (irow, csf_count, "IROW", "initialize_sparse_eigenvectors")
        read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
        !      IF (myid .EQ. 0) PRINT *, 'iniestsd ...........'
        if (csf_count /= file_csf_count .or.  mpi_rank /= file_mpi_rank                      &
            .or. file_mpi_process_count /= mpi_process_count)                        &
            stop "initialize_sparse_eigenvectors: restart file CSF count mismatch"

        !          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
        !          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.

        do j = mpi_rank + 1, ns, mpi_process_count
            joff = (j*(j-1))/2
            read (restart_matrix_unit) nelc,elsto,(hmx(ir),ir=1,nelc),                  &
                (irow(ir),ir=1,nelc)
            hmx(nelc) = hmx(nelc) - hamiltonian_energy_shift ! Shift the diagonal
            do ir = 1, nelc
                ap(irow(ir) + joff) = hmx(ir)
            end do
        end do

        ! Let each node have a complete copy of ap

        call allreduce_mpi_real_sum (ap, (ns*(ns+1))/2)

        ! To be in step with other cases, go through the whole block.
        !
        ! This is not necessary since currently the file pointer is moved
        ! to the absolute position of the .res files which is always counted
        ! from the begining of the .res files of each node. Besides, the
        ! following segment seems not working properly for the last block.
        ! Xinghong He 98-12-14

        !mylast = j - nprocs
        !DO j = mylast, ncf, nprocs
        !   READ (restart_matrix_unit)
        !ENDDO

        call dalloc (hmx,"HMX", "initialize_sparse_eigenvectors")
        call dalloc (irow, "irow", "initialize_sparse_eigenvectors")

        call alloc (eigval,ns,"EIGVAL","initialize_sparse_eigenvectors")
        call alloc (vec,ns*initial_vector_count,"VEC","initialize_sparse_eigenvectors")
        call alloc (real_workspace,8*ns,"real_workspace","initialize_sparse_eigenvectors")
        call alloc (integer_workspace,5*ns,"integer_workspace","initialize_sparse_eigenvectors" )
        call alloc (ifail,ns, "IFAIL","initialize_sparse_eigenvectors" )

        !      CALL DSPEVX ('Vectors also','In a range','Upper triangular',     &
        !                NS,AP,-1.0D0,-1.0D0,1,NIV,0.d0,                        &
        !                NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
        abstol = 2*slamch("S")
        call dspevx ("V","I","U",                                     &
            ns,ap,-1.0_real64,-1.0_real64,1,initial_vector_count,abstol,                         &
            nfound,eigval,vec,ns,real_workspace,integer_workspace,ifail,info)
        ierr = -abs (info)

        !******************************************************************

        !       ..Build the Basis.

        call initialize_real_vector (csf_count*initial_vector_count, 0.0_real64, basis, 1)
        !       ...scatter the vectors
        do j = 1, initial_vector_count
            call dcopy (ns, vec(ns*(j-1)+1:ns*(j-1)+(ns)),1, basis(csf_count*(j-1)+1), 1)
        end do
        call dcopy (initial_vector_count, eigval,1,basis(initial_vector_count*csf_count+1),1)

        call dalloc (ap, "AP", "initialize_sparse_eigenvectors")
        call dalloc (eigval,"EIGVAL", "initialize_sparse_eigenvectors")
        call dalloc (vec, "VEC", "initialize_sparse_eigenvectors")
        call dalloc (real_workspace, "real_workspace", "initialize_sparse_eigenvectors")
        call dalloc (integer_workspace, "integer_workspace", "initialize_sparse_eigenvectors")
        call dalloc (ifail, "IFAIL", "initialize_sparse_eigenvectors")

        return
    end subroutine initialize_sparse_eigenvectors
end module initialize_sparse_eigenvectors_interface
