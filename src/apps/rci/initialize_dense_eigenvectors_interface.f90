module initialize_dense_eigenvectors_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine initialize_dense_eigenvectors (nmax, csf_count, initial_vector_count, basis, hmx)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in):: nmax, csf_count, initial_vector_count
            real(real64) :: basis(*), hmx(*)
        end subroutine  initialize_dense_eigenvectors
    end interface
end module initialize_dense_eigenvectors_interface
