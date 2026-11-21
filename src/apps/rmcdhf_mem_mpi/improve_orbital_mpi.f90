module improve_orbital_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: improve_orbital_mpi
contains

    !***********************************************************************
    subroutine improve_orbital_mpi (extended_optimal_level, orbital_index, sort_by_self_consistency, maximum_orbital_damping)
        ! Module dependencies
        !   The difference from the serial version is that it calls MPI        *
        !   version subroutines (setlagmpi, cofpotmpi, matrixmpi, newcompi).   *
        !                                                                      *
        !   Improve the orbital J.                                             *
        !                                                                      *
        !   Call(s) to: [RSCF92]: DACON, DAMPCK, DAMPOR, LAGCON, matrixmpi,    *
        !                         newcompi, ORTHOR, ROTATE, SETCOF, setlagmpi, *
        !                         SOLVE, XPOT, YPOT.                           *
        !               [LIB92]: ORTHSC, QUAD.                                 *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
        !   Modified by Xinghong He                 Last update: 05 Aug 1988   *
        !   Modified for ifort -i8 by A. Kramida (AK) Last update 22 Mar 2016  *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc, ralloc
        use scf_damping_factors, only: orbital_damping_factors
        use atomic_calculation_state, only: working_tolerance, orbital_updates_per_scf_iteration
        use radial_grid_state, only: radial_point_count, radial_grid_derivative
        use radial_dirac_workspace, only: trial_orbital_last_grid_index, trial_large_component, trial_large_origin_parameter, &
            trial_small_component
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_gather, mpi_integer, &
            mpi_process_count, mpi_rank
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_origin_exponents, orbital_principal_n, &
            orbital_symmetry_labels
        use scf_orthogonality_control, only: enforce_orbital_orthogonality
        use scf_iteration_state, only: radial_solution_methods, orbital_self_consistency_errors, direct_coefficient_count, &
            direct_coefficient_capacity, direct_coefficients, direct_orbital_indices
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters

        ! Procedure interfaces
        use build_orbital_potential_mpi_interface, only: build_orbital_potential_mpi
        use calculate_deferred_corrections_interface, only: calculate_deferred_corrections
        use solve_radial_orbital_interface, only: solve_radial_orbital

        ! Module dependencies
        use radial_orbital_orthogonalization, only: orthogonalize_radial_orbitals

        ! Procedure interfaces
        !GG      USE matrixmpi_I
        !GG      USE newcompi_I
        use calculate_lagrange_multipliers_mpi_interface, only: calculate_lagrange_multipliers_mpi

        ! Module dependencies
        use radial_quadrature, only: integrate_radial_grid

        ! Procedure interfaces
        use evaluate_radial_self_consistency_interface, only: evaluate_radial_self_consistency
        use check_orbital_damping_interface, only: check_orbital_damping
        use damp_orbital_update_interface, only: damp_orbital_update
        use orthogonalize_same_kappa_orbitals_interface, only: orthogonalize_same_kappa_orbitals
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        real(real64), intent(inout) :: maximum_orbital_damping
        logical, intent(in) :: extended_optimal_level
        logical :: sort_by_self_consistency

        ! Local constants
        real(real64), parameter :: method_switch_norm_error_threshold = 2.0e-1_real64
        real(real64), parameter :: method_switch_energy_change_threshold = 5.0e-3_real64
        real(real64), parameter :: iteration_reduction_energy_threshold = 1.0e-4_real64

        ! Local declarations
        integer :: ipr, potential_grid_point_count, orbital_sign_reversal_count, matching_grid_index, radial_node_count, i, nwww
        integer :: i_mpi, ndcof_max, ntot, i_last, iproc, ndcip, jproc
        integer :: ifound, ind_buf, k
        real(real64) :: previous_relative_energy_change, origin_exponent, old_orbital_energy, radial_norm_squared, &
            normalization_factor, relative_energy_change, normalization_error, orbital_damping_factor
        logical :: solve_failed, first
        real(real64), pointer :: da_buffer(:)
        integer, pointer :: nda_buffer(:),ndcof_buffer(:)

        ! Initial values
        !-----------------------------------------------
        !
        !   C Froese Fischer's IPR and ED1 parameter
        !
        data ipr/ 0/
        data previous_relative_energy_change/ 0.0_real64/
        data first/ .false./
        !
        !
        !-----------------------------------------------------------------------
        !AK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
        !      ISIZE = sizeof(NCF)
        i_mpi = mpi_integer
        !      if (ISIZE.EQ.8) I_MPI = MPI_INTEGER8
        !
        origin_exponent = orbital_origin_exponents(orbital_index)
        !
        !   C Froese Fischer's parameters IPR, ED1, ED2 are set and
        !   used in this routine and in DAMPCK
        !
        do
            old_orbital_energy = orbital_energy_parameters(orbital_index)
            !
            !   Set up the exchange potential and arrays XU, XV as appropriate
            !
            !   Set coefficients for YPOT, XPOT, DACON
            !   Compute direct potential, exchange potential
            !   Add in Lagrange-multiplier contribution
            !   Add in derivative-terms contribution
            !
            potential_grid_point_count = radial_point_count
            call build_orbital_potential_mpi (extended_optimal_level, orbital_index, potential_grid_point_count)
            !
            !   Calculate deferred corrections
            !
            call calculate_deferred_corrections (orbital_index)
            !
            !   Solve the Dirac equation
            !
            !AK     NDCOF, NDA, and DA cannot be reduced by finding the max of NDCOF
            !          and summing up DA from all nodes. The correct replacement is below
            !cjb MPI_ALLREDUCE NDCOF
            !      call MPI_ALLREDUCE(ndcof,ndcof_buffer,1,
            !     :          I_MPI,MPI_MAX,MPI_COMM_WORLD,ierror)
            !      if (ndcof_buffer .gt. ndcof) then
            !        print *, 'improvmpi: ndcof, ndcof_buffer, myid=',
            !     &                     ndcof, ndcof_buffer, myid
            !        write (*,'(A35,4I10)') 'improv11: myid,NDCOF,ndcbfr=',
            !     &    myid,NDCOF,ndcof_buffer
            !        do indcof = ndcof+1, ndcof_buffer
            !           da(indcof) = 0.0
            !           nda(indcof) = 0
            !        enddo
            !        ndcof = ndcof_buffer
            !      endif
            orbital_sign_reversal_count = 0
            !cjb MPI_ALLREDUCE DA
            !      if (ndcof.gt.0) then
            !         call alloc(pda_buffer, ndcof, 8)
            !         call alloc(pnda_buffer, ndcof, ISIZE)
            !         call MPI_ALLREDUCE(da,da_buffer,ndcof,
            !     :          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
            !         da(1:ndcof) = da_buffer(1:ndcof)
            !      nda(1:ndcof) = nda_buffer(1:ndcof)
            !         call dalloc(pda_buffer)
            !         call dalloc(pnda_buffer)
            !      ENDIF
            !AK      Consolidate NDA and DA from all nodes
            if (mpi_process_count>1) then
                call alloc(ndcof_buffer,mpi_process_count,"ndcof_buffer","improve_orbital_mpi")
                call mpi_gather(direct_coefficient_count, 1, i_mpi, ndcof_buffer,            &
                    1, i_mpi, 0, mpi_comm_world, mpi_error_code)
                call mpi_bcast (ndcof_buffer, mpi_process_count, i_mpi, 0,           &
                    mpi_comm_world, mpi_error_code)
                ndcof_max = 0
                do i = 1, mpi_process_count
                    if (ndcof_buffer(i)>ndcof_max) ndcof_max=ndcof_buffer(i)
                end do
                if (ndcof_max>0) then
                    ntot = ndcof_max*mpi_process_count
                    call alloc(nda_buffer,ntot,"nda_buffer","improve_orbital_mpi")
                    call mpi_gather(direct_orbital_indices, ndcof_max, i_mpi, nda_buffer,     &
                        ndcof_max, i_mpi, 0, mpi_comm_world, mpi_error_code)
                    call mpi_bcast (nda_buffer, ntot, i_mpi, 0,            &
                        mpi_comm_world, mpi_error_code)
                    call alloc(da_buffer, ntot, "da_buffer","improve_orbital_mpi")
                    call mpi_gather(direct_coefficients, ndcof_max, mpi_double_precision,   &
                        da_buffer, ndcof_max, mpi_double_precision, 0,       &
                        mpi_comm_world, mpi_error_code)
                    call mpi_bcast (da_buffer, ntot, mpi_double_precision, &
                        0, mpi_comm_world, mpi_error_code)
                    i_last = 0
                    do iproc = 1, mpi_process_count
                        ndcip = ndcof_buffer(iproc)
                        do jproc = 1, ndcip
                            ifound = 0
                            ind_buf = (iproc - 1)*ndcof_max + jproc
                            do k = 1, i_last
                                if (nda_buffer(ind_buf) == direct_orbital_indices(k)) then
                                    direct_coefficients(k) = direct_coefficients(k) + da_buffer(ind_buf)
                                    ifound = k
                                    exit
                                end if
                            end do
                            if (ifound==0) then
                                i_last = i_last + 1
                                do while (i_last > direct_coefficient_capacity)
                                    if (direct_coefficient_capacity > 0) then
                                        direct_coefficient_capacity = 2*direct_coefficient_capacity
                                        !GG                           CALL alcsca (PNTNDA, PNTRDA, NDDIM, 2)
                                        call ralloc(direct_orbital_indices,direct_coefficient_capacity,"direct_orbital_indices", &
                                            "improve_orbital_mpi")
                                    else
                                        direct_coefficient_capacity = 64
                                        !GG                           CALL alcsca (PNTNDA, PNTRDA, NDDIM, 1)
                                        call alloc(direct_orbital_indices,direct_coefficient_capacity,"direct_orbital_indices", &
                                            "improve_orbital_mpi")
                                    end if
                                end do
                                direct_orbital_indices(i_last) = nda_buffer(ind_buf)
                                direct_coefficients(i_last) = da_buffer(ind_buf)
                            end if
                        end do
                    end do
                    direct_coefficient_count = ndcof_max
                    call dalloc(nda_buffer,"nda_buffer","improve_orbital_mpi")
                    call dalloc(da_buffer,"da_buffer","improve_orbital_mpi")
                end if
                call dalloc(ndcof_buffer,"ndcof_buffer","improve_orbital_mpi")
            end if

            call solve_radial_orbital (orbital_index, solve_failed, orbital_sign_reversal_count, matching_grid_index, &
                & radial_node_count)
            !
            !   Upon failure issue message; take corrective action if possible
            !
            if (solve_failed) then
                if (mpi_rank == 0) write (*, '(/," Failure; equation for orbital ",1i2,1a2,&
                    &" could not be solved using method ",1i1)') &
                    orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                        radial_solution_methods(orbital_index)
                if (radial_solution_methods(orbital_index) /= 2) then
                    radial_solution_methods(orbital_index) = 2
                    !XHH orthsc does not have any argument
                    !    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
                    call orthogonalize_radial_orbitals
                    !CFF        ... avoid rediagonalization
                    !            IF (EOL) THEN
                    !               CALL MATRIXMPI
                    !               CALL NEWCOMPI (WTAEV)
                    !            ENDIF
                    call calculate_lagrange_multipliers_mpi (extended_optimal_level)
                    cycle
                else
                    if (mpi_rank == 0) write (*, '(/,/," ****** Error in SUBROUTINE improve_orbital_mpi ******"/,&
                        &" Convergence not obtained"/)')
                    !CALL TIMER (0)
                    stop
                end if
            end if
            !
            !   Compute norm of radial function
            !
            radial_workspace_a(1) = 0.0_real64
            radial_workspace_a(2:trial_orbital_last_grid_index) = (trial_large_component(2:trial_orbital_last_grid_index)**2+ &
                trial_small_component(2:trial_orbital_last_grid_index)**2)*radial_grid_derivative(2:trial_orbital_last_grid_index)
            integration_last_index = trial_orbital_last_grid_index

            call integrate_radial_grid (radial_norm_squared)

            !   Determine self-consistency [multiplied by SQRT(UCF(J))]

            call evaluate_radial_self_consistency (orbital_index)
            !
            !   Normalize
            !
            normalization_factor = 1.0_real64/dsqrt(radial_norm_squared)
            trial_large_origin_parameter = trial_large_origin_parameter*normalization_factor
            trial_large_component(:trial_orbital_last_grid_index) = trial_large_component(:trial_orbital_last_grid_index)* &
                normalization_factor
            trial_small_component(:trial_orbital_last_grid_index) = trial_small_component(:trial_orbital_last_grid_index)* &
                normalization_factor
            !
            !   Check if different method should be used or if improvement
            !   count should be reduced
            !
            relative_energy_change = dabs(1.0_real64 - old_orbital_energy/orbital_energy_parameters(orbital_index))
            if (radial_solution_methods(orbital_index) == 1) then
                normalization_error = dmax1(dabs(1.0_real64 - dsqrt(radial_norm_squared)),dabs(normalization_factor - 1.0_real64))
                if (relative_energy_change<method_switch_energy_change_threshold .and. &
                    & normalization_error>method_switch_norm_error_threshold) then
                    radial_solution_methods(orbital_index) = 2
                    cycle
                end if
            else
                if (relative_energy_change<iteration_reduction_energy_threshold .and. orbital_updates_per_scf_iteration>1) &
                    orbital_updates_per_scf_iteration = orbital_updates_per_scf_iteration - 1
            end if
            exit
        end do
        !
        !   Damp the orbital --- if not converged
        !
        if (orbital_self_consistency_errors(orbital_index) > working_tolerance) then
            call check_orbital_damping (ipr, orbital_index, previous_relative_energy_change, old_orbital_energy)
            orbital_damping_factor = dabs(orbital_damping_factors(orbital_index))
        else
            orbital_damping_factor = 0.0_real64                           ! take the whole new orbital
        end if
        call damp_orbital_update (orbital_index, orbital_sign_reversal_count, orbital_damping_factor)

        !   Orthogonalize all orbitals of the same kappa in the order
        !   fixed, spectroscopic, correlation orbitals. The order of
        !   orbitals in the latter two classes are sorted according
        !   to their self-consistency and energy.

        if (enforce_orbital_orthogonality) then
            !CALL orthor (J, inv)
            nwww = orbital_count
            call orthogonalize_same_kappa_orbitals (nwww, orbital_index, sort_by_self_consistency)
        end if
        !
        !   Print details of iteration
        !
        if (mpi_rank == 0)                                                &
            write (*, '(1x,1i2,1a2,1p,1d16.7,1x,1i2,d11.3,1d10.2,1d10.2,0p,f6.3,1x,1i5,1x,1i5,1x,1i2,1x,1i2)') &
            orbital_principal_n(orbital_index),orbital_symmetry_labels(orbital_index),orbital_energy_parameters(orbital_index), &
            & radial_solution_methods(orbital_index), &
            & orbital_origin_parameters(orbital_index),orbital_self_consistency_errors(orbital_index), &
            !cjb DNORM-1 -> SQRT(DNORM)-1
            !cjb                    DNORM - 1, ODAMPJ, JP, MF(J), INV, NNP
            sqrt(radial_norm_squared)-1, orbital_damping_factor, matching_grid_index, orbital_last_grid_index(orbital_index), &
            & orbital_sign_reversal_count, radial_node_count
        maximum_orbital_damping = dmax1(maximum_orbital_damping,dabs(orbital_damping_factor))

        return
    end subroutine improve_orbital_mpi
end module improve_orbital_mpi_interface
