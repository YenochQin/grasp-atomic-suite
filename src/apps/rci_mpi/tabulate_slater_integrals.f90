!***********************************************************************
!                                                                      *
subroutine tabulate_slater_integrals (mpi_rank, mpi_process_count, n, max_orbital_two_j)
    ! Module dependencies
    !                                                                      *
    !  Input:
    !     myid, nprocs
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
    use distributed_slater_integrals, only: indteirk, valteirk, rk_local_count
    use coulomb_rank_offsets, only: kmax, kstart
    use orbital_csf_basis, only: orbital_count, orbital_l, orbital_two_j
    use angular_symmetry_labels, only: left_justified_j_labels => jlbl, parity_labels => jlbp, right_justified_j_labels => jlbr
    use slater_radial_integral, only: slater

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
    kstart(0) = 1
    !
    !   Find 2*JMAX; it should not be greater than PARAMETER KMAX
    !
    max_orbital_two_j = orbital_two_j(1)
    do i = 2, orbital_count
        if (orbital_two_j(i) > max_orbital_two_j) max_orbital_two_j = orbital_two_j(i)
    end do

    if (max_orbital_two_j > kmax) then
        stop "genintrk: KMAX too small"
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
                                        indteirk(n) = ((orbital_index*key+ib)*key+ic)*key+id
                                        valteirk(n) = slater (orbital_index,ib,ic,id,k)
                                    end if
                                end if
                            end do
                        end if
                    end do
                end do
            end do
            kstart(k+1) = n + 1
        end do
        !
        !     Allocate memory for integral book keeping
        !
        if (gen) exit
        call alloc (indteirk,n,"INDTEIRK", "GENINTRK")
        call alloc (valteirk,n,"VALTEIRK", "GENINTRK")

        ! Initialization is necessary in the mpi version

        do i = 1, n
            indteirk(i) = 0
            valteirk(i) = 0.0_real64
        end do

        if (mpi_rank == 0)                                          &
            print *, "Allocating space for ",n," Rk integrals"

        gen = .true.
    end do

    return
end subroutine tabulate_slater_integrals
