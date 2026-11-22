! Form starting eigenvector estimates using a dense Hamiltonian submatrix.
module initialize_dense_eigenvectors_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_dense_eigenvectors
contains

    !     Routine for providing initial estimates from the diagonal
    !     of the matrix. This way was used by Dvdson in atomic structure
    !     calculations. It should be used to obtain estimates when nothing
    !     else is available.
    !
    ! Block version by Xinghong He. Last revision: 18 Jun 1998
    subroutine initialize_dense_eigenvectors (nmax, csf_count, initial_vector_count, basis, hmx)

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc

        ! Procedure interfaces
        use blas_lapack_interfaces, only: dcopy, dspevx

        ! Module dependencies
        use real_vector_initialization, only: initialize_real_vector
        implicit none (type, external)

        ! Arguments
        integer, intent(in):: nmax, csf_count, initial_vector_count
        real(real64) :: basis(*), hmx(*)

        ! Local declarations
        integer :: mpi_rank, mpi_process_count, info, ierr, ns, i, nfound, j
        real(real64), pointer :: ap(:), eigval(:), vec(:), real_workspace(:)
        integer, pointer :: integer_workspace(:), ifail(:)

        mpi_rank = 0
        mpi_process_count = 1

        ns = min (nmax, csf_count)

        call alloc (ap, ns*(ns+1)/2, "AP", "initialize_dense_eigenvectors")
        call initialize_real_vector (ns*(ns+1)/2, 0.0_real64, ap, 1)

        ! Get the upper left sub-matrix

        do i = 1, (ns*(ns+1))/2
            ap(i) = hmx(i)
        end do

        call alloc (eigval,   ns, "EIGVAL", "initialize_dense_eigenvectors")
        call alloc (vec,   ns*initial_vector_count, "VEC", "initialize_dense_eigenvectors")
        call alloc (real_workspace,  8*ns,  "real_workspace", "initialize_dense_eigenvectors")
        call alloc (integer_workspace, 5*ns,  "integer_workspace", "initialize_dense_eigenvectors")
        call alloc (ifail,     ns,  "IFAIL", "initialize_dense_eigenvectors")

        call dspevx ("Vectors also","In a range","Upper triangular",      &
            &          ns,ap,-1.0_real64,-1.0_real64,1,initial_vector_count,0.0_real64,                               &
            &          nfound,eigval,vec,ns,real_workspace,integer_workspace,ifail,info)
        ierr = -abs (info)

        ! Build the Basis.

        call initialize_real_vector (csf_count*initial_vector_count, 0.0_real64, basis, 1)

        ! scatter the vectors

        do j = 1, initial_vector_count
            call dcopy (ns, vec(ns*(j-1)+1:ns*(j-1)+(ns)),1, basis(csf_count*(j-1)+1), 1)
        end do

        call dcopy (initial_vector_count, eigval, 1, basis(initial_vector_count*csf_count+1), 1)

        call dalloc (ap, "AP", "initialize_dense_eigenvectors")
        call dalloc (eigval, "EIGVAL", "initialize_dense_eigenvectors")
        call dalloc (vec, "VEC", "initialize_dense_eigenvectors")
        call dalloc (real_workspace, "real_workspace", "initialize_dense_eigenvectors")
        call dalloc (integer_workspace, "integer_workspace", "initialize_dense_eigenvectors")
        call dalloc (ifail, "IFAIL", "initialize_dense_eigenvectors")

        return
    end subroutine initialize_dense_eigenvectors
end module initialize_dense_eigenvectors_interface
