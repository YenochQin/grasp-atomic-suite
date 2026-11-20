module hamiltonian_correction_flags
    implicit none (type, external)

    ! Module state
    logical :: use_first_order_csf_perturbation, include_transverse_interaction, include_vacuum_polarization, &
        estimate_self_energy, include_normal_mass_shift, include_specific_mass_shift
end module hamiltonian_correction_flags
