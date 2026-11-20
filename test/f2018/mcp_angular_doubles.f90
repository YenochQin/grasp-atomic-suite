! Predictable coefficients isolate the actual MCP reader and potential assembly.
module csf_density_matrix_element_interface
    implicit none (type, external)
    private
    public :: csf_density_matrix_element
contains
    function csf_density_matrix_element(extended_optimal_level, bra_csf_index, ket_csf_index, symmetry_block_index) result(value)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)
        logical, intent(in) :: extended_optimal_level
        integer, intent(in) :: bra_csf_index, ket_csf_index
        integer :: symmetry_block_index
        real(real64) :: value
        value = 1.0_real64
    end function csf_density_matrix_element
end module csf_density_matrix_element_interface

module direct_angular_coefficient_interface
    implicit none (type, external)
    private
    public :: direct_angular_coefficient
contains
    function direct_angular_coefficient(multipole_rank, csf_index, orbital_a, orbital_b) result(value)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)
        integer, intent(in) :: multipole_rank, orbital_a, orbital_b
        integer :: csf_index
        real(real64) :: value
        value = 0.0_real64
    end function direct_angular_coefficient
end module direct_angular_coefficient_interface

module exchange_angular_coefficient_interface
    implicit none (type, external)
    private
    public :: exchange_angular_coefficient
contains
    function exchange_angular_coefficient(multipole_rank, csf_index, orbital_a, orbital_b) result(value)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)
        integer :: multipole_rank, csf_index
        integer, intent(in) :: orbital_a, orbital_b
        real(real64) :: value
        value = 0.0_real64
    end function exchange_angular_coefficient
end module exchange_angular_coefficient_interface
