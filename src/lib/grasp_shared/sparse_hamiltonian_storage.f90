module sparse_hamiltonian_storage
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64, int64
    implicit none (type, external)

    ! Module interface
    private :: real64, int64

    ! Module state
    ! Each stored column segment ends with its diagonal. End offsets are cumulative within
    ! the local arrays. MPI ranks store only their assigned columns.
    integer(int64) :: stored_matrix_element_count
    real(real64), pointer :: hamiltonian_values(:)
    integer, pointer :: stored_segment_end_offsets(:)
    integer, pointer :: matrix_partner_indices(:)
    integer :: ntpitmp(6)
    integer(int64) :: nelmnttmp
    integer :: ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp,  &
        ncoretmp, nvpitmp, nkeitmp, nvintitmp, ncftmp
    real(real64) :: cutofftmp
end module sparse_hamiltonian_storage
