module load_final_transition_basis_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_final_transition_basis
contains

    !***********************************************************************
    !                                                                      *
    subroutine load_final_transition_basis(final_core_subshell_count, state_name)
        ! Module dependencies
        !                                                                      *
        !   Open, check, load data from and close the  .csl  file. This file   *
        !   is always attached to stream 21.                                   *
        !                                                                      *
        !   Call(s) to: [LIB92]: IQ, LENGTH, LODCSL, OPENFL.                   *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use atomic_calculation_state, only: electron_count, final_electron_count
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_kappa, orbital_l, orbital_principal_n, &
            orbital_symmetry_labels, orbital_two_j
        use biorthogonal_orbital_sets, only: final_orbital_symmetry_labels, final_csf_count, final_orbital_count, &
            final_orbital_principal_n, final_orbital_kappa, final_orbital_l, final_orbital_two_j

        use csf_full_basis_reader, only: load_csf_basis
        use transition_basis_blocks, only: transition_block_count, cumulative_csf_counts, final_transition_block_count, &
            final_cumulative_csf_counts
        implicit none (type, external)

        ! Arguments
        integer , intent(out) :: final_core_subshell_count
        character , intent(in) :: state_name*24

        ! Local declarations
        integer :: io_status, ncorer
        !-----------------------------------------------
        !
        !
        !     Common relevant for the radial wave functions of the final state.
        !
        !
        !   The  .csl  file is FORMATTED; it must exist
        !
        call load_csf_basis(trim(state_name),ncorer,io_status)
        if (io_status /= 0) error stop "Transition descriptor load failed"
        if (transition_block_count > size(final_cumulative_csf_counts)) error stop "Transition block capacity exceeded"
        final_transition_block_count=transition_block_count
        final_cumulative_csf_counts=0
        final_cumulative_csf_counts(:transition_block_count)=cumulative_csf_counts(:transition_block_count)
        !
        !   Check if the core should be redefined
        !
        final_core_subshell_count = ncorer
        !      DO 3 I = NCORER+1,NW
        !         IFULLI = NKJ(I)+1
        !         DO 2 J = 1,NCF
        !    2    CONTINUE
        !         CALL CONVRT (NP(I),RECORD,LENTH)
        !         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
        !         PRINT *, ' in all CSFs; including this'
        !         PRINT *, ' subshell in the core;'
        !         NCORE = NCORE+1
        !    3 CONTINUE
        !
        final_electron_count = electron_count
        final_orbital_count = orbital_count
        final_csf_count = csf_count
        final_orbital_symmetry_labels(:orbital_count) = orbital_symmetry_labels(:orbital_count)
        final_orbital_principal_n(:orbital_count) = orbital_principal_n(:orbital_count)
        final_orbital_kappa(:orbital_count) = orbital_kappa(:orbital_count)
        final_orbital_l(:orbital_count) = orbital_l(:orbital_count)
        final_orbital_two_j(:orbital_count) = orbital_two_j(:orbital_count)
        return
    end subroutine load_final_transition_basis
end module load_final_transition_basis_interface
