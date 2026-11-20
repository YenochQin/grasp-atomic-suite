module basis_csf_two_j_plus_one_interface
    implicit none (type, external)
    private
    public :: basis_csf_two_j_plus_one
contains
integer function basis_csf_two_j_plus_one(index)
    implicit none (type, external)
    integer,intent(in) :: index
    if(index<1.or.index>3) error stop "Unexpected level index"
    basis_csf_two_j_plus_one=2*index-1
end function basis_csf_two_j_plus_one
end module basis_csf_two_j_plus_one_interface
