! Access the intermediate coupling of the final-state CSF in the active basis.
module right_csf_intermediate_coupling_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: right_csf_intermediate_coupling
contains

    integer function right_csf_intermediate_coupling(loc, csf_index)
        use csf_quantum_queries, only: left_intermediate_coupling
        integer, intent(in) :: loc, csf_index

        right_csf_intermediate_coupling = left_intermediate_coupling(loc,csf_index)
    end function right_csf_intermediate_coupling
end module right_csf_intermediate_coupling_interface
