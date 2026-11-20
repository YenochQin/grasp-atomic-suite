module basis_csf_two_j_plus_one_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: basis_csf_two_j_plus_one
contains

    integer function basis_csf_two_j_plus_one(csf_index)
        ! Module dependencies
        use rmcdhf_csf_input, only: rmcdhf_csf_symmetry
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: csf_index
        basis_csf_two_j_plus_one = abs(rmcdhf_csf_symmetry(csf_index))
    end function basis_csf_two_j_plus_one
end module basis_csf_two_j_plus_one_interface
