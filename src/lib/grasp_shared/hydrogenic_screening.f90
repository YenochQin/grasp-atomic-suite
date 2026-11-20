module hydrogenic_screening
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: orbital_screening_constants(max_orbitals)
end module hydrogenic_screening
