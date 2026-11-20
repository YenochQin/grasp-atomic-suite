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
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use isotope_properties,        only: sqn, dmomnm, qmomb,                  &
        hfsi=>sqn, hfsd=>dmomnm, hfsq=>qmomb
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0

    ! Procedure interfaces
    use write_hyperfine_energies_interface, only: write_hyperfine_energies
    implicit none (type, external)
    !-----------------------------------------------
    !
    !    Write nuclear data
    !
    write (summary_output_unit, '("Nuclear spin                        ",1p,d22.15," au")') hfsi
    write (summary_output_unit, '("Nuclear magnetic dipole moment      ",1p,d22.15," n.m.")') hfsd
    write (summary_output_unit, '("Nuclear electric quadrupole moment  ",1p,d22.15," barns")') hfsq
    write (hyperfine_output_unit, '("Nuclear spin                        ",1p,d22.15," au")') hfsi
    write (hyperfine_output_unit, '("Nuclear magnetic dipole moment      ",1p,d22.15," n.m.")') hfsd
    write (hyperfine_output_unit, '("Nuclear electric quadrupole moment  ",1p,d22.15," barns")') hfsq
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
