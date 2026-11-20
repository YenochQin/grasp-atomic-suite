module block_mean_energies
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64), pointer :: mean_energy_by_block(:)
end module block_mean_energies
