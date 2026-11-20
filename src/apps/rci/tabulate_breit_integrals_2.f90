module tabulate_breit_integrals_2_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: tabulate_breit_integrals_2
contains

    !***********************************************************************
    !                                                                      *
    subroutine tabulate_breit_integrals_2 (mpi_rank, mpi_process_count, nb, max_orbital_two_j)
        ! Module dependencies
        !                                                                      *
        !  Input:
        !     MPI rank and process count
        !  Output:
        !   N - Number of integrals
        !   j2max - max of 2*J
        !
        !       Generate the list of Breit Integrals of type 2                 *
        !       that could arise from a set of orbitals                        *
        !       of orbitals.                                                   *
        !                                                                      *
        !       This routine is similar to genintrk                            *
        !                                                                      *
        !     Written by Per Jonsson            October 2014                   *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use array_memory_management, only: alloc
        use breit_integral_tables, only: breit_integral_indices_2, breit_integral_values_2
        use orbital_csf_basis, only: orbital_two_j, orbital_count
        use breit_rank_offsets, only: maximum_breit_multipole_rank, second_breit_integral_group_offsets

        ! Procedure interfaces
        use lookup_breit_radial_integral_interface, only: lookup_breit_radial_integral
        use satisfies_breit_group2_triangle_interface, only: satisfies_breit_group2_triangle

        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: mpi_rank, mpi_process_count
        integer             :: nb, max_orbital_two_j

        ! Local declarations
        logical :: gen
        integer :: i, key, k, orbital_index, ib, ic, id
        !-----------------------------------------------
        !

        key = orbital_count + 1
        !   Type 2 angular coefficients use rank nu+1 (nu >= 0).
        !   The empty rank-zero interval and rank one start at the same offset.
        second_breit_integral_group_offsets(0:1) = 1
        !
        !   Find 2*JMAX; it should not be greater than PARAMETER KMAX
        !
        max_orbital_two_j = orbital_two_j(1)
        do i = 2, orbital_count
            if (orbital_two_j(i) > max_orbital_two_j) max_orbital_two_j = orbital_two_j(i)
        end do

        if (max_orbital_two_j > maximum_breit_multipole_rank) then
            stop "genintbreit2: KMAX too small"
        end if
        !
        !   Make the breit integrals. Note that there are no symmetry relations that
        !   can be utilized. See paper by Grasnt.
        !   When GEN is false, sweep through to find dimension
        !
        gen = .false.

        do
            nb = 0
            do k = 1, max_orbital_two_j
                do orbital_index = 1, orbital_count
                    do ib = 1, orbital_count
                        do ic = 1, orbital_count
                            do id = 1, orbital_count
                                if (satisfies_breit_group2_triangle(orbital_index,ib,ic,id,k)) then
                                    nb = nb + 1
                                    if (gen .and.(mod(nb,mpi_process_count) == mpi_rank)) then
                                        breit_integral_indices_2(nb) = ((orbital_index*key+ib)*key+ic)*key+id
                                        breit_integral_values_2(nb) = lookup_breit_radial_integral (2,orbital_index,ib,ic,id,k)
                                    end if
                                end if
                            end do
                        end do
                    end do
                end do
                second_breit_integral_group_offsets(k+1) = nb + 1
            end do
            !
            !     Allocate memory for integral book keeping
            !
            if (gen) exit
            call alloc (breit_integral_indices_2,nb,"INDT2","genintbreit2")
            call alloc (breit_integral_values_2,nb,"VALT2","genintbreit2")

            ! Initialization is necessary in the mpi version

            do i = 1, nb
                breit_integral_indices_2(i) = 0
                breit_integral_values_2(i) = 0.0_real64
            end do

            if (mpi_rank == 0)                                        &
                print *, "Computing",nb," Breit integrals of type 2"

            gen = .true.
        end do

        return
    end subroutine tabulate_breit_integrals_2
end module tabulate_breit_integrals_2_interface
