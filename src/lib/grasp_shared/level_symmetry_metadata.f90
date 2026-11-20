! Stores twice-J and parity for the active, final, and initial sets of levels.
module level_symmetry_metadata
    implicit none (type, external)

    ! Module constants
    integer, parameter::nblk0 = 50  ! Capacity of the symmetry-block lookup array.

    ! Module state
    integer, pointer :: iatjpo(:), iaspar(:)
    integer, pointer :: iatjpoff(:), iasparff(:)
    integer, pointer :: iatjpoii(:), iasparii(:)
    integer :: jpgg(nblk0)
end module level_symmetry_metadata
