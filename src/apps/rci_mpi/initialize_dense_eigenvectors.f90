module initialize_dense_eigenvectors_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_dense_eigenvectors
contains

    !************************************************************************
    !
    subroutine initialize_dense_eigenvectors (nmax, csf_count, initial_vector_count, basis, hmx)

        ! Module dependencies
        !     Routine for providing initial estimates from the diagonal
        !     of the matrix. This way was used by Dvdson in atomic structure
        !     calculations. It should be used to obtain estimates when nothing
        !     else is available.
        !
        !   Block version by Xinghong He            Last revision: 18 Jun 1998
        !
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc
        use mpi_state, only: mpi_error_code, mpi_process_count, mpi_rank

        ! Procedure interfaces
        use blas_lapack_interfaces, only: dcopy, dspevx

        ! Module dependencies
        use real_vector_initialization, only: initialize_real_vector

        ! Procedure interfaces
        use blas_lapack_interfaces, only: slamch
        implicit none (type, external)

        ! Arguments
        integer, intent(in):: nmax, csf_count, initial_vector_count
        real(real64) :: basis(*), hmx(*)

        ! Local declarations
        integer :: info, ns, i, nfound, j, ihmx, iap
        real(real64), pointer :: ap(:), eigval(:), vec(:), real_workspace(:)
        integer, pointer :: integer_workspace(:), ifail(:)
        real(real64) :: abstol
        !-----------------------------------------------------------------------
        ns = min (nmax, csf_count)

        call alloc (ap, ns*(ns+1)/2, "AP", "INIESTDM")
        call initialize_real_vector (ns*(ns+1)/2, 0.0_real64, ap, 1)

        !  Get the upper left sub-matrix
        !  ap is global, hmx is local !!!

        ihmx = 0
        do i = mpi_rank + 1, ns, mpi_process_count
            iap = i*(i-1)/2
            do j = 1, i
                ap(j+iap) = hmx(j+ihmx)
            end do
            ihmx = ihmx + i
        end do

        call alloc (eigval,   ns, "EIGVAL", "INIESTDM")
        call alloc (vec,   ns*initial_vector_count, "VEC", "INIESTDM")
        call alloc (real_workspace,  8*ns,  "WORK", "INIESTDM")
        call alloc (integer_workspace, 5*ns,  "IWORK", "INIESTDM")
        call alloc (ifail,     ns,  "IFAIL", "INIESTDM")

        !      CALL DSPEVX ('Vectors also','In a range','Upper triangular',     &
        !     &          NS,AP,-1.0D0,-1.0D0,1,NIV,0.d0,                        &
        !     &          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
        abstol = 2*slamch("S")
        call dspevx ("V","I","U",                                        &
            ns,ap,-1.0_real64,-1.0_real64,1,initial_vector_count,abstol,                            &
            nfound,eigval,vec,ns,real_workspace,integer_workspace,ifail,info)
        mpi_error_code = -abs (info)

        !  Build the Basis.

        call initialize_real_vector (csf_count*initial_vector_count, 0.0_real64, basis, 1)

        !  scatter the vectors

        do j = 1, initial_vector_count
            call dcopy (ns, vec(ns*(j-1)+1:ns*(j-1)+(ns)),1, basis(csf_count*(j-1)+1), 1)
        end do

        call dcopy (initial_vector_count, eigval, 1, basis(initial_vector_count*csf_count+1), 1)

        call dalloc (ap, "AP", "INIESTDM")
        call dalloc (eigval, "EIGVAL", "INIESTDM")
        call dalloc (vec, "VEC", "INIESTDM")
        call dalloc (real_workspace, "WORK", "INIESTDM")
        call dalloc (integer_workspace, "IWORK", "INIESTDM")
        call dalloc (ifail, "IFAIL", "INIESTDM")

        return
    end subroutine initialize_dense_eigenvectors
end module initialize_dense_eigenvectors_interface
