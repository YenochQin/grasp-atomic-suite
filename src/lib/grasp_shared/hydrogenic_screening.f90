module hydrogenic_screening
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: sigma(max_orbitals)
end module hydrogenic_screening
