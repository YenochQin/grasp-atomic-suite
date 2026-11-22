! Check whether the two CSF occupations admit a one-electron matrix element.
module one_electron_occupation_allowed_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: one_electron_occupation_allowed
contains

    !   This routine is to check the occupation condition for one electron *
    !   operator.
    !
    ! Yu Zou. Last revision: 8/16/00
    integer function one_electron_occupation_allowed (ja, jb)
        ! Module dependencies
        use orbital_csf_basis, only: orbital_count
        use csf_quantum_queries, only: csf_occupation
        implicit none (type, external)

        ! Arguments
        integer  :: ja
        integer  :: jb

        ! Local declarations
        integer :: k, orbital_index, iqa, iqb

        one_electron_occupation_allowed = 0
        k = 0
        do orbital_index = 1, orbital_count
            iqa = csf_occupation(orbital_index,ja)
            iqb = csf_occupation(orbital_index,jb)
            if (iqa == iqb) cycle
            k = k + 1
            if (k > 2) return
            if (iabs(iqa - iqb) <= 1) cycle
            return
        end do
        if (k==2 .or. k==0) one_electron_occupation_allowed = 1
        return
    end function one_electron_occupation_allowed
end module one_electron_occupation_allowed_interface
