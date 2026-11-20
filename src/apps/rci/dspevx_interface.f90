module dspevx_interface
    ! Procedure interfaces
    ! Re-export the canonical LP64 LAPACK interface.
    use blas_lapack_interfaces, only: dspevx
    implicit none (type, external)
end module dspevx_interface
