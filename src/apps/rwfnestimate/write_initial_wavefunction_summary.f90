module write_initial_wavefunction_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_initial_wavefunction_summary
contains

    !    Prints  a summary of the complete  list of subshell radial wave-
    !    functions on NUNIT.
    !
    ! Written by Farid A. Parpia. Last revision: 15 Dec 1992
    !
    ! Modified by Yanting Li, output <r> for each orbital, 19 Aug 2022
    subroutine write_initial_wavefunction_summary(nunit)
        !  Module dependencies
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_origin_exponents, orbital_principal_n, &
            orbital_symmetry_labels
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters
        use initial_orbital_source, only: orbital_source_description
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: nunit

        ! Local declarations
        integer :: orbital_index, lenth

        write (nunit, '("Shell",6x,"e",11x,"p0",8x,"gamma",8x,"<r>",6x,"MTP","  SRC"/)')
        do orbital_index = 1, orbital_count
            lenth = len_trim(orbital_source_description(orbital_index))
            ! --------Output P(2) and Q(2) for each orbital
            ! --------Output <r> (RINT(I,I,1)) for each orbital
            write (nunit, '(1x,i2,a2,4d12.4,i5,2x,a3)') &
                orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), orbital_energy_parameters( &
                & orbital_index), &
                & orbital_origin_parameters(orbital_index), orbital_origin_exponents(orbital_index), &
                    & integrate_weighted_radial_overlap(orbital_index, &
                & orbital_index,1)&
                , orbital_last_grid_index(orbital_index), orbital_source_description(orbital_index)(1:lenth)
        end do
        return
        ! -Format for P(2) and Q(2)
        ! -Format for <r>
        return

    end subroutine write_initial_wavefunction_summary
end module write_initial_wavefunction_summary_interface
