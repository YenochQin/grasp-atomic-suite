module tabulate_slater_integrals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: tabulate_slater_integrals
contains

    !***********************************************************************
    !                                                                      *
    subroutine tabulate_slater_integrals (mpi_rank, mpi_process_count, n, max_orbital_two_j)
        ! Module dependencies
        !                                                                      *
        !  Input:
        !     MPI rank and process count
        !  Output:
        !   N - Number of integrals
        !   j2max - max of 2*J
        !
        !       Generate the list of Rk Integrals that could arise from a set  *
        !       of orbitals.                                                   *
        !                                                                      *
        !     Written by Per Jonsson                                           *
        !   MPI version by Xinghong He            Last revision: 22 Jun 1998   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc
        use distributed_slater_integrals, only: local_slater_integral_indices, local_slater_integral_values
        use coulomb_rank_offsets, only: maximum_coulomb_multipole_rank, coulomb_integral_offsets_by_rank
        use orbital_csf_basis, only: orbital_count, orbital_l, orbital_two_j
        use slater_radial_integral, only: calculate_slater_radial_integral

        ! Procedure interfaces
        use satisfies_coulomb_triangle_interface, only: satisfies_coulomb_triangle
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: mpi_rank, mpi_process_count
        integer, intent(out) :: n, max_orbital_two_j

        ! Local declarations
        logical :: gen
        integer :: key, i, k, orbital_index, ib, ic, id
        !-----------------------------------------------------------------------

        key = orbital_count + 1
        coulomb_integral_offsets_by_rank(0) = 1
        !
        !   Find 2*JMAX; it should not be greater than PARAMETER KMAX
        !
        max_orbital_two_j = orbital_two_j(1)
        do i = 2, orbital_count
            if (orbital_two_j(i) > max_orbital_two_j) max_orbital_two_j = orbital_two_j(i)
        end do

        if (max_orbital_two_j > maximum_coulomb_multipole_rank) then
            stop "tabulate_slater_integrals: maximum multipole rank is too small"
        end if
        !
        !   Make the RK integrals: IA <= IB, IA <= IC, IA <= ID, IB <= ID
        !   When GEN is false, sweep through to find dimension
        !
        gen = .false.

        do
            n = 0
            do k = 0, max_orbital_two_j
                do orbital_index = 1, orbital_count
                    do ib = orbital_index, orbital_count
                        do ic = orbital_index, orbital_count
                            if (satisfies_coulomb_triangle(orbital_l(orbital_index),k,orbital_l(ic))) then
                                do id = ib, orbital_count
                                    if (satisfies_coulomb_triangle(orbital_l(ib),k,orbital_l(id))) then
                                        n = n + 1
                                        if (gen .and. (mod(n,mpi_process_count) == mpi_rank)) then
                                            local_slater_integral_indices(n) = ((orbital_index*key+ib)*key+ic)*key+id
                                            local_slater_integral_values(n) = calculate_slater_radial_integral &
                                                (orbital_index,ib,ic,id,k)
                                        end if
                                    end if
                                end do
                            end if
                        end do
                    end do
                end do
                coulomb_integral_offsets_by_rank(k+1) = n + 1
            end do
            !
            !     Allocate memory for integral book keeping
            !
            if (gen) exit
            call alloc (local_slater_integral_indices,n,"local_slater_integral_indices", "tabulate_slater_integrals")
            call alloc (local_slater_integral_values,n,"local_slater_integral_values", "tabulate_slater_integrals")

            ! Initialization is necessary in the mpi version

            do i = 1, n
                local_slater_integral_indices(i) = 0
                local_slater_integral_values(i) = 0.0_real64
            end do

            if (mpi_rank == 0)                                               &
                &      print *, "Allocating space for ",n," Rk integrals"

            gen = .true.
        end do

        return
    end subroutine tabulate_slater_integrals
end module tabulate_slater_integrals_interface
