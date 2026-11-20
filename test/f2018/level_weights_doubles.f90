module basis_csf_two_j_plus_one_interface
    implicit none (type, external)
    private
    public :: basis_csf_two_j_plus_one
contains
    integer function basis_csf_two_j_plus_one(csf_index)
        implicit none (type, external)
        integer, intent(in) :: csf_index
        if (csf_index < 1 .or. csf_index > 3) error stop "Unexpected level index"
        basis_csf_two_j_plus_one = 2*csf_index - 1
    end function basis_csf_two_j_plus_one
end module basis_csf_two_j_plus_one_interface
