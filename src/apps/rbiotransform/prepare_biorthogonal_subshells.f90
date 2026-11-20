module prepare_biorthogonal_subshells_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: prepare_biorthogonal_subshells
contains

    !***********************************************************************
    !                                                                      *
    subroutine prepare_biorthogonal_subshells(ntestg, ncore1, ncore2)
        ! Module dependencies
        !                                                                      *
        !   This subroutine determines the number of kappa quantum numbers     *
        !   KAMAX together with the number of orbitals of each kappa.          *
        !                                                                      *
        !   This subroutine also checks if the orbital ordering is normal      *
        !   or reversed. If reversed then JA and JB have to be permuted        *
        !   In the subroutine ti1tv                                            *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !***********************************************************************
        use grasp_parameters,   only: max_orbitals
        use biorthogonal_orbital_sets,         only: initial_orbital_count, initial_orbital_kappa, &
            initial_orbital_principal_n, final_orbital_count, final_orbital_kappa, final_orbital_principal_n
        use biorthogonal_subshell_layout,         only: max_biorthogonal_kappa_groups, initial_core_counts_by_group, &
            final_core_counts_by_group, initial_orbital_counts_by_group, final_orbital_counts_by_group, kappa_by_group, &
            initial_kappa_group_by_orbital, final_kappa_group_by_orbital, initial_orbitals_by_group, final_orbitals_by_group, &
            initial_orbital_positions_by_group, final_orbital_positions_by_group, kappa_group_count
        use biorthogonal_orbital_order, only: final_orbital_order_code, initial_orbital_order_code
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: ntestg
        integer, intent(in) :: ncore1
        integer, intent(in) :: ncore2

        ! Local declarations
        integer :: isort(2*max_orbitals)
        integer :: ntestl, ntest, i, j, k, ital, nref
        logical :: klar
        !-----------------------------------------------
        !
        !
        !
        ntestl = 0
        ntest = max0(ntestl,ntestg)

        isort(:2*max_orbitals) = 100

        kappa_by_group(:max_biorthogonal_kappa_groups) = 0
        initial_orbital_counts_by_group(:max_biorthogonal_kappa_groups) = 0
        final_orbital_counts_by_group(:max_biorthogonal_kappa_groups) = 0
        initial_core_counts_by_group(:max_biorthogonal_kappa_groups) = 0
        final_core_counts_by_group(:max_biorthogonal_kappa_groups) = 0

        initial_kappa_group_by_orbital(:max_orbitals) = 0
        final_kappa_group_by_orbital(:max_orbitals) = 0

        initial_orbitals_by_group(:max_biorthogonal_kappa_groups,:max_biorthogonal_kappa_groups) = 0
        final_orbitals_by_group(:max_biorthogonal_kappa_groups,:max_biorthogonal_kappa_groups) = 0

        initial_orbital_positions_by_group(:max_biorthogonal_kappa_groups,:max_orbitals) = 0
        final_orbital_positions_by_group(:max_biorthogonal_kappa_groups,:max_orbitals) = 0
        !
        !   Sort the kappa quantum numbers
        !
        do k = 1, initial_orbital_count + final_orbital_count
            if (k <= initial_orbital_count) then
                ital = initial_orbital_kappa(k)
            else
                ital = final_orbital_kappa(k-initial_orbital_count)
            end if
            i = k - 1
            klar = .false.
            do while (i>0 .and. .not.klar)
                if (ital <= isort(i)) then
                    isort(i+1) = isort(i)
                    i = i - 1
                else
                    klar = .true.
                end if
            end do
            isort(i+1) = ital
        end do
        !
        !   Determine the unique set of kappa IKAPPA
        !
        kappa_group_count = 1
        kappa_by_group(1) = isort(1)
        do k = 1, 2*max_orbitals - 1
            if (isort(k)==isort(k+1) .or. isort(k+1)>=100) cycle
            kappa_group_count = kappa_group_count + 1
            kappa_by_group(kappa_group_count) = isort(k+1)
        end do
        !
        !  Make a connection between each kappa and a number in the
        !  range [1,KAMAX] as to know on which file to dump the data
        !  Determine the number of shells NSHLII for each I in the
        !  range [1,KAMAX]
        !
        if (ntest >= 10) then
            write (*, *) "******************"
            write (*, *) " Entering kapdata"
            write (*, *) "******************"
            write (*, *)
            write (*, *) "There are", initial_orbital_count, "orbitals in the initial state"
            write (*, *) "with the following n and kappa quantum numbers"
        end if

        do j = 1, initial_orbital_count
            if (ntest >= 10) write (*, *) "orbital number", j, "n and kappa", initial_orbital_principal_n&
                (j), initial_orbital_kappa(j)
            if (j <= ncore1) then
                do i = 1, kappa_group_count
                    if (kappa_by_group(i) /= initial_orbital_kappa(j)) cycle
                    initial_kappa_group_by_orbital(j) = i
                    initial_core_counts_by_group(i) = initial_core_counts_by_group(i) + 1
                    initial_orbital_counts_by_group(i) = initial_orbital_counts_by_group(i) + 1
                    initial_orbitals_by_group(i,initial_orbital_counts_by_group(i)) = j
                    initial_orbital_positions_by_group(i,j) = initial_orbital_counts_by_group(i)
                end do
            else
                do i = 1, kappa_group_count
                    if (kappa_by_group(i) /= initial_orbital_kappa(j)) cycle
                    initial_kappa_group_by_orbital(j) = i
                    initial_orbital_counts_by_group(i) = initial_orbital_counts_by_group(i) + 1
                    initial_orbitals_by_group(i,initial_orbital_counts_by_group(i)) = j
                    initial_orbital_positions_by_group(i,j) = initial_orbital_counts_by_group(i)
                end do
            end if
        end do

        if (ntest >= 10) then
            write (*, *) "There are", final_orbital_count, "orbitals in the final state"
            write (*, *) "with the following n and kappa quantum numbers"
        end if
        do j = 1, final_orbital_count
            if (ntest >= 10) write (*, *) "orbital number", j, "n and kappa=", &
                final_orbital_principal_n(j), final_orbital_kappa(j)
            if (j <= ncore2) then
                do i = 1, kappa_group_count
                    if (kappa_by_group(i) /= final_orbital_kappa(j)) cycle
                    final_kappa_group_by_orbital(j) = i
                    final_core_counts_by_group(i) = final_core_counts_by_group(i) + 1
                    final_orbital_counts_by_group(i) = final_orbital_counts_by_group(i) + 1
                    final_orbitals_by_group(i,final_orbital_counts_by_group(i)) = j
                    final_orbital_positions_by_group(i,j) = final_orbital_counts_by_group(i)
                end do
            else
                do i = 1, kappa_group_count
                    if (kappa_by_group(i) /= final_orbital_kappa(j)) cycle
                    final_kappa_group_by_orbital(j) = i
                    final_orbital_counts_by_group(i) = final_orbital_counts_by_group(i) + 1
                    final_orbitals_by_group(i,final_orbital_counts_by_group(i)) = j
                    final_orbital_positions_by_group(i,j) = final_orbital_counts_by_group(i)
                end do
            end if
        end do

        if (ntest >= 10) then
            write (*, *) "Total number of different kappa", kappa_group_count
            do i = 1, kappa_group_count
                write (*, *) "L=", i, "corresponds to kappa=", kappa_by_group(i)
                write (*, *) "nr of init.  orb. with this kappa=", initial_orbital_counts_by_group(i)
                write (*, *) "nr of final. orb. with this kappa=", final_orbital_counts_by_group(i)
            end do
            do i = 1, kappa_group_count
                write (*, *) "Position in initial state list"
                do j = 1, initial_orbital_counts_by_group(i)
                    write (*, *) "L=", i, "orb. nr", j, ",position", initial_orbitals_by_group(i,j)
                end do
                write (*, *) "Position in final   state list"
                do j = 1, final_orbital_counts_by_group(i)
                    write (*, *) "L=", i, "orb. nr", j, ",position", final_orbitals_by_group(i,j)
                end do
                write (*, *) "Relative positions for initial state orbitals"
                do j = 1, initial_orbital_count
                    if (initial_orbital_positions_by_group(i,j) == 0) cycle
                    write (*, *) "Orbital", j, "is nr", initial_orbital_positions_by_group(i,j), "with kappa"&
                        , kappa_by_group(i)
                end do
                write (*, *) "Relative positions for final state orbitals"
                do j = 1, final_orbital_count
                    if (final_orbital_positions_by_group(i,j) == 0) cycle
                    write (*, *) "Orbital", j, "is nr", final_orbital_positions_by_group(i,j), "with kappa"&
                        , kappa_by_group(i)
                end do
            end do

        end if
        !
        !  Check if the orbital ordering is normal or reversed.
        !
        !  For normal ordering
        !
        !    NPII(NSHLPII(I,1)) < NPII(NSHLPII(I,2)) < ... < NPII(NSHLPII(I,NSHLII(I))
        !
        !  For reversed ordering
        !
        !    NPII(NSHLPII(I,1)) > NPII(NSHLPII(I,2)) > ... > NPII(NSHLPII(I,NSHLII(I))
        !
        initial_orbital_order_code = 0
        do i = 1, kappa_group_count
            nref = 0
            do j = 1 + initial_core_counts_by_group(i), initial_orbital_counts_by_group(i)
                if (initial_orbital_principal_n(initial_orbitals_by_group(i,j)) < nref) initial_orbital_order_code = 1
                nref = initial_orbital_principal_n(initial_orbitals_by_group(i,j))
            end do
        end do

        final_orbital_order_code = 0
        do i = 1, kappa_group_count
            nref = 0
            do j = 1 + final_core_counts_by_group(i), final_orbital_counts_by_group(i)
                if (final_orbital_principal_n(final_orbitals_by_group(i,j)) < nref) final_orbital_order_code = 1
                nref = final_orbital_principal_n(final_orbitals_by_group(i,j))
            end do
        end do

        if (initial_orbital_order_code /= final_orbital_order_code) then
            write (*, *) " Orbital order of the initial and final states"
            write (*, *) " should be the same. STOP"
            stop
        end if
        !
        !   If not normal order check if reversed order
        !
        if (initial_orbital_order_code == 1) then
            do i = 1, kappa_group_count
                nref = 0
                do j = initial_orbital_counts_by_group(i), 1 + initial_core_counts_by_group(i), -1
                    if (initial_orbital_principal_n(initial_orbitals_by_group(i,j)) < nref) initial_orbital_order_code = 2
                    nref = initial_orbital_principal_n(initial_orbitals_by_group(i,j))
                end do
            end do

            do i = 1, kappa_group_count
                nref = 0
                do j = final_orbital_counts_by_group(i), 1 + final_core_counts_by_group(i), -1
                    if (final_orbital_principal_n(final_orbitals_by_group(i,j)) < nref) final_orbital_order_code = 2
                    nref = final_orbital_principal_n(final_orbitals_by_group(i,j))
                end do
            end do
        end if

        if (initial_orbital_order_code==2 .or. final_orbital_order_code==2) then
            write (*, *) " The orbital order is neither normal or reversed"
            write (*, *) " STOP"
            stop
        end if

        !w      write(*,*) 'Give nordii'
        !w      read(*,*) nordii
        !w      nordff = nordii

        if (initial_orbital_order_code == 0) then
            write (*, *) " Normal orbital ordering"
        else
            write (*, *) " Reverse orbital ordering"
        end if

        write (*, *) "*****************"
        write (*, *) " Leaving kapdata"
        write (*, *) "*****************"
        write (*, *)

        return
    end subroutine prepare_biorthogonal_subshells
end module prepare_biorthogonal_subshells_interface
