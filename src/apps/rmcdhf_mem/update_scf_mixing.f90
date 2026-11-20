!***********************************************************************
!                                                                      *
subroutine update_scf_mixing(sum)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This routine computes the level weights, the generalized occupa-   *
    !   tion numbers,  and average energy for  EOL calculations;  this     *
    !   information and the eigenvectors are then printed out.             *
    !                                                                      *
    !   Call(s) to: [LIB92]: IQ.                                           *
    !               [RSCF92]: CSFWGT, DSUBRS.                              *
    !                                         Last revision: 24 Dec 1992   *
    !   Block version by Xinghong He          Last revision: 05 Aug 1998   *
    !   Midified by G. Gaigalas                              05 Feb 2017   *
    !      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
    !                                  JCUPA(NNNW*NCF)                     *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use atomic_calculation_state, only: iccmin, ncmin, weight, wt
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use orbital_csf_basis, only: max_orbitals, orbital_count
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use level_symmetry_metadata,          only: jpgg
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
    use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
        csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
    use block_mean_energies, only: eavblk
    use transition_block_indices, only: idxblk
    use serial_mpi_context, only: mpi_process_count, mpi_rank

    ! Procedure interfaces
    use csf_density_matrix_element_interface, only: csf_density_matrix_element

    ! Module dependencies
    use csf_quantum_queries, only: iq

    ! Procedure interfaces
    use write_csf_weights_interface, only: write_csf_weights
    implicit none (type, external)

    ! Arguments
    real(real64), intent(out) :: sum

    ! Local declarations
    integer :: j, symmetry_block_index, csf_index, noff, jall, jblock
    real(real64) :: weitj, ee
    logical :: eol
    !-----------------------------------------------
    !
    !     POINTER (pncfblk, ncfblk(0:*))
    !
    !
    !   Compute weighting factors
    !
    sum = 0.0_real64
    do j = 1, ncmin
        weitj = weight(j)
        if (weitj == (-2.0_real64)) then
            wt(j) = 1.0_real64
        else if (weitj == (-1.0_real64)) then
            !GGGG
            !GG            WT(J) = IATJPO(J)
            jblock = idxblk(j)     ! Block number of this state
            wt(j) = abs(jpgg(jblock))
        else
            wt(j) = weight(j)
        end if
        sum = sum + wt(j)
    end do

    wt(:ncmin) = wt(:ncmin)/sum
    !
    !   Compute generalised occupation numbers
    !   <----- Distributed ----->
    !
    eol = .true.

    do j = 1, orbital_count
        sum = 0.0_real64
        do symmetry_block_index = 1, symmetry_block_count
            do csf_index = mpi_rank + 1, csf_counts_by_block(symmetry_block_index), mpi_process_count
                sum = sum + csf_density_matrix_element(eol,csf_index,csf_index,symmetry_block_index)*iq(j, &
                    & csf_index + csf_offsets_by_block(symmetry_block_index))
            end do
        end do
        ucf(j) = sum
    end do
    !
    !   Write out level energies and weights
    !
    write (*, '(/,"Optimise on the following level(s):"/)')
    sum = 0.0_real64
    noff = 0
    do jall = 1, ncmin
        symmetry_block_index = idxblk(jall)                       ! Block number of this state
        ee = eavblk(symmetry_block_index) + shifted_eigenvalues(jall)
        write (*, '("Level ",1i2,4x,"Energy = ",1p,1d19.12,4x,"Weight = ",1d12.5)') iccmin(jall), ee, wt(jall)
        if (ldbpg(5)) then
            call ensure_debug_output()
            write (angular_debug_unit, *) jall, symmetry_block_index, csf_counts_by_block(symmetry_block_index), &
                & mixing_offsets_by_block(symmetry_block_index)
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"Configuration mixing coefficients:")')
            call ensure_debug_output()
            write (angular_debug_unit, '(1x,1p,6d12.4)') (csf_mixing_coefficients(csf_index + noff),csf_index=1, &
                & csf_counts_by_block(symmetry_block_index))
            noff = noff + csf_counts_by_block(symmetry_block_index)
        end if
        sum = sum + wt(jall)*ee
    end do

    call write_csf_weights (.true.)
    !
    !   Write out average energy
    !
    if (ncmin > 1) write (*, '(/,"Weighted average energy of these levels = ",1p,d18.10)') sum
    !
    !   Write out generalized occupation numbers
    !
    write (*, '(/,"Generalised occupation numbers:"/)')
    write (*, '(1x,1p,6d12.4)') (ucf(csf_index),csf_index=1,orbital_count)


    return
end subroutine update_scf_mixing
