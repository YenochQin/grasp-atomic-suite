! Return the parity of a CSF in the loaded MCDHF basis.
module basis_csf_parity_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: basis_csf_parity
contains

    integer function basis_csf_parity(csf_index)
        ! Module dependencies
        use rmcdhf_csf_input, only: rmcdhf_csf_symmetry
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: csf_index
        basis_csf_parity = sign(1,rmcdhf_csf_symmetry(csf_index))
    end function basis_csf_parity
end module basis_csf_parity_interface
