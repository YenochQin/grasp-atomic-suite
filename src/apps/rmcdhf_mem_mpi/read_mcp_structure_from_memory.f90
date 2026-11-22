module read_mcp_structure_from_memory_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mcp_structure_from_memory
contains

    ! Written by  G. Gaigalas               Vilnius, September 2021
    subroutine read_mcp_structure_from_memory(jblock)
        ! Module dependencies
        use sparse_hamiltonian_storage,     only: matrix_partner_indices
        use rmcdhf_mcp_memory, only: imax2_30, imin2_30, irow_30
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: jblock

        ! Local declarations
        integer :: i

        do i = imin2_30(jblock)+1, imax2_30(jblock)
            matrix_partner_indices(i-imin2_30(jblock)) = irow_30(i)
        end do

        return
    end subroutine read_mcp_structure_from_memory
end module read_mcp_structure_from_memory_interface
