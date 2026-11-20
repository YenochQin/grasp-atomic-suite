module write_input_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_input_summary
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_input_summary
        ! Module dependencies
        use calculation_file_units, only: hyperfine_output_unit, summary_output_unit
        !                                                                      *
        !   Generates the first part of  hfs92.sum  (on stream 24 and 29).     *
        !                                                                      *
        !   Call(s) to: ENGOUTH                                                *
        !                                                                      *
        !   Written by P. Jonsson                 Last revision: 20 Oct 1999   *
        !                                                                      *
        !***********************************************************************
        use isotope_properties,        only: nuclear_spin, nuclear_magnetic_dipole_moment, nuclear_electric_quadrupole_moment

        ! Procedure interfaces
        use write_hyperfine_energies_interface, only: write_hyperfine_energies
        implicit none (type, external)
        !-----------------------------------------------
        !
        !    Write nuclear data
        !
        write (summary_output_unit, '("Nuclear spin                        ",1p,d22.15," au")') nuclear_spin
        write (summary_output_unit, '("Nuclear magnetic dipole moment      ",1p,d22.15," n.m.")') nuclear_magnetic_dipole_moment
        write (summary_output_unit, '("Nuclear electric quadrupole moment  ",1p,d22.15," barns")') &
            nuclear_electric_quadrupole_moment
        write (hyperfine_output_unit, '("Nuclear spin                        ",1p,d22.15," au")') nuclear_spin
        write (hyperfine_output_unit, '("Nuclear magnetic dipole moment      ",1p,d22.15," n.m.")') nuclear_magnetic_dipole_moment
        write (hyperfine_output_unit, '("Nuclear electric quadrupole moment  ",1p,d22.15," barns")') &
            nuclear_electric_quadrupole_moment
        !
        !    Write the list of eigenpair indices
        !
        !      CALL ENGOUTH (EAV, EVAL, IATJPO, IASPAR, IVEC, NVEC, 3)
        !
        return
        !
        return
        !
    end subroutine write_input_summary
end module write_input_summary_interface
