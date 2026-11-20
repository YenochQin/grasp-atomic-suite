module write_initial_wavefunction_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_initial_wavefunction_summary
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_initial_wavefunction_summary(nunit)
        ! Module dependencies
        !                                                                      *
        !   Prints  a summary of the complete  list of subshell radial wave-   *
        !   functions on NUNIT.                                                *
        !                                                                      *
        !   Call(s) to: [LIB92]:                                               *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        !
        !...Modified by Yanting Li, output <r> for each orbital, 19 Aug 2022
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
        !-----------------------------------------------
        !
        !
        write (nunit, '("Shell",6x,"e",11x,"p0",8x,"gamma",8x,"<r>",6x,"MTP","  SRC"/)')
        !
        do orbital_index = 1, orbital_count
            lenth = len_trim(orbital_source_description(orbital_index))
            !--------Output P(2) and Q(2) for each orbital
            !         WRITE (NUNIT, 301) NP(I), NH(I), E(I), PZ(I), GAMA(I), PF(2,I), QF(2,I&
            !            ), MF(I), SOURCE(I)(1:LENTH)
            !--------Output <r> (RINT(I,I,1)) for each orbital
            write (nunit, '(1x,i2,a2,4d12.4,i5,2x,a3)') &
                orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), orbital_energy_parameters( &
                & orbital_index), &
                & orbital_origin_parameters(orbital_index), orbital_origin_exponents(orbital_index), &
                    & integrate_weighted_radial_overlap(orbital_index, &
                & orbital_index,1)&
                , orbital_last_grid_index(orbital_index), orbital_source_description(orbital_index)(1:lenth)
            !         WRITE (NUNIT,302) SOURCE(I)(1:LENTH)
        end do
        !
        return
        !
        !-Format for P(2) and Q(2)
        !  300 FORMAT('Shell',6X,'e',11X,'p0',8X,'gamma',8X,'P(2)',7X,'Q(2)',6X,'MTP',&
        !         '  SRC'/)
        !  301 FORMAT(1X,I2,A2,5D12.4,I5,2X,A3)
        !-Format for <r>
        return
        !  302 FORMAT ('      Source: ',A)
        !
    end subroutine write_initial_wavefunction_summary
end module write_initial_wavefunction_summary_interface
