module write_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_biorthogonal_mcp
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_biorthogonal_mcp(name, startup_directory, ik, ntestg, inpci)
        ! Module dependencies
        use biorthogonal_subshell_layout, only: kappa_coefficient_units
        use biorthogonal_subshell_layout, only: biorthogonal_angular_unit
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   This routine controls the computation  and storage of the values   *
        !   and all indices of the angular coefficients                        *
        !                                                                      *
        !                                                                      *
        !                   T  (ab)                                            *
        !                    rs                                                *
        !                                                                      *
        !   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
        !   are orbital sequence numbers.  r and s are configuration           *
        !   state function indices.                                            *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, DALLOC, TNSRJJ                         *
        !               [GENMCP]: QSORT.                                       *
        !                                                                      *
        !   Written by Per Jonsson                Last revision: JUne 1996     *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_orbitals, orbital_index_radix
        use array_memory_management, only: dalloc
        !      USE def_C
        !      USE EIGV_C
        use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes
        use transition_basis_blocks, only: cumulative_csf_counts, transition_block_count
        use orbital_csf_basis,           only: csf_count, orbital_count, orbital_kappa, csf_subshell_occupations
        use biorthogonal_transition_selection, only: itjqj1, itjqj2, jqj1, jqj2, ntrans
        use biorthogonal_orbital_order, only: final_orbital_order_code, initial_orbital_order_code
        use biorthogonal_subshell_layout,         only: kappa_group_count, max_biorthogonal_kappa_groups, &
            initial_kappa_group_by_orbital, final_kappa_group_by_orbital
        use mpi_state, only: mpi_process_count, mpi_rank
        use one_electron_scalar_coefficients, only: calculate_one_electron_scalar_coefficients
        use closed_shell_mcp, only: calculate_closed_shell_mcp
        use csf_quantum_queries, only: csf_total_j_code

        ! Procedure interfaces
        use open_biorthogonal_angular_cache_interface, only: open_biorthogonal_angular_cache
        use sort_biorthogonal_mcp_interface, only: sort_biorthogonal_mcp
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: ik
        integer, intent(in) :: ntestg
        integer  :: inpci
        character, intent(in) :: name*24
        character  :: startup_directory*128

        ! Local constants
        integer, parameter :: nvmax = 100
        integer, parameter :: orbital_label_radix_squared = orbital_index_radix*orbital_index_radix
        real(real64), parameter :: cutoff = 1.0e-10_real64

        ! Local declarations
        character(len=len(name)) :: cache_state_name
        integer :: llistt(max_biorthogonal_kappa_groups)
        character(len=24) :: coefficient_filename
        integer :: nakinv(max_orbitals)
        integer :: i, ntestl, ntest, ka, iopar, j, ncf0, iblk, k, &
            bra_csf_index, ket_csf_index, ibb, orbital_index, ket_orbital_index, lab, jan, jbn, l
        real(real64) :: shell_coefficients(max_orbitals)
        logical   :: comp, avail
        !-----------------------------------------------
        !
        ! Locals
        !     POINTER (pscr, scr(1))
        !     POINTER (pciout, ciout(1))
        !

        write (output_unit, *) "NBLOCK,(NCFBLK(i),i=1,NBLOCK)"
        write (output_unit, *) transition_block_count, (cumulative_csf_counts(i),i=1,transition_block_count)

        ntestl = 0
        ntest = max(ntestg,ntestl)
        ntest = 0
        !
        !   Set the rank (zero) and parity (even) for the one-particle
        !   coefficients
        !
        ka = 0
        iopar = 1
        !
        !   Check if angular data is available. If available read this data.
        !   If not available calculate the data
        !
        ! The cache helper accepts a writable name and only reads its value here.
        cache_state_name = name
        call open_biorthogonal_angular_cache(cache_state_name, avail, kappa_group_count)
        write (output_unit, *) "AVAIL=", avail
        if (avail) return

        if (ik == 1) then
            nakinv(:max_orbitals) = initial_kappa_group_by_orbital(:max_orbitals)
        else
            nakinv(:max_orbitals) = final_kappa_group_by_orbital(:max_orbitals)
        end if



        write (output_unit, *) " open sorted ang. file .TB(NF)", 200
        j = index(name," ")
        open(newunit=biorthogonal_angular_unit, file=name(1:j-1)//".TB", status="UNKNOWN", form=&
            "UNFORMATTED", position="asis")

        ncf0 = 1
        do iblk = 1, transition_block_count
            !
            !   Open scratchfiles to dump the T coefficients for each kappa
            !
            do k = 1, kappa_group_count
                write(coefficient_filename, "(A,I0)") "fort.", 80 + k
                open(newunit=kappa_coefficient_units(k), file=coefficient_filename, &
                    status="UNKNOWN", form="UNFORMATTED", position="asis")
            end do
            !
            !   Initialize the counters for the total number of T coefficients
            !
            llistt(:max_biorthogonal_kappa_groups) = 0
            !
            !   JA and JB respectively refer to the initial and final states
            !   in the list of NCF configurations
            !
            do bra_csf_index = mpi_rank+ncf0,cumulative_csf_counts(iblk),mpi_process_count
                if (mod(bra_csf_index,100)==0 .and. ik==1) then
                    write (*, *) " JA1 =", bra_csf_index, bra_csf_index - ncf0 + 1
                else if (mod(bra_csf_index,100)==0 .and. ik==2) then
                    write (*, *) " JA2 =", bra_csf_index, bra_csf_index - ncf0 + 1
                end if
                !
                do ket_csf_index = ncf0, cumulative_csf_counts(iblk)
                    !
                    !   Call the MCT package to compute T coefficients
                    !
                    if (ntrans == 1) then
                        comp = .false.
                        if (ik == 1) then
                            do ibb = 1, jqj1
                                if (csf_total_j_code(bra_csf_index) /= itjqj1(ibb)) cycle
                                comp = .true.
                            end do
                        else
                            do ibb = 1, jqj2
                                if (csf_total_j_code(bra_csf_index) /= itjqj2(ibb)) cycle
                                comp = .true.
                            end do
                        end if
                    else
                        comp = .true.
                    end if


                    if (.not.comp) cycle
                    !            write(*,*) JA,JB
                    !GG               CALL TNSRJJ (KA, IOPAR, JA, JB, IA, IB, TSHELL)
                    call calculate_one_electron_scalar_coefficients(bra_csf_index,ket_csf_index,orbital_index,ket_orbital_index, &
                        & shell_coefficients)
                    if (orbital_index == 0) cycle
                    if (orbital_index == ket_orbital_index) then
                        do orbital_index = 1, orbital_count
                            !
                            !   If T coefficient is greater than zero and the kappa quantum numbers
                            !   of the two orbitals are the same dump to file
                            !   In a later version use a buffer with a reasonable record length
                            !
                            if (dabs(shell_coefficients(orbital_index)) <= cutoff) cycle
                            llistt(nakinv(orbital_index)) = llistt(nakinv(orbital_index)) + 1
                            lab = orbital_index*orbital_index_radix + orbital_index
                            write (kappa_coefficient_units(nakinv(orbital_index))) bra_csf_index - ncf0 + 1, ket_csf_index - ncf0 &
                                & + 1, lab&
                                , shell_coefficients(orbital_index)
                        end do
                    else
                        if (dabs(shell_coefficients(1))>cutoff .and. orbital_kappa(orbital_index)==orbital_kappa( &
                            & ket_orbital_index)) then
                            llistt(nakinv(orbital_index)) = llistt(nakinv(orbital_index)) + 1
                            if (initial_orbital_order_code==0 .and. final_orbital_order_code==0) then
                                !
                                !   Experssion for normal orbital ordering
                                !
                                lab = orbital_index*orbital_index_radix + ket_orbital_index
                                jan = bra_csf_index - ncf0 + 1
                                jbn = ket_csf_index - ncf0 + 1
                            else if (initial_orbital_order_code==1 .and. final_orbital_order_code==1) then
                                !
                                !   Experssion for reversed orbital ordering
                                !
                                lab = ket_orbital_index*orbital_index_radix + orbital_index
                                jan = ket_csf_index - ncf0 + 1
                                jbn = bra_csf_index - ncf0 + 1
                            else
                                write (*, *) "SOMETHING WRONG"
                                stop
                            end if
                            write (kappa_coefficient_units(nakinv(orbital_index))) jan, jbn, lab, shell_coefficients(1)
                        end if
                    end if
                    !
                end do
            end do


            !
            !   sort the MCP data into inegral based lists.
            !
            do l = 1, kappa_group_count
                if (llistt(l) > 0) then
                    call sort_biorthogonal_mcp (l, llistt(l), ik, name, kappa_group_count)
                else
                    if (l == 1) write (biorthogonal_angular_unit) csf_count, orbital_count, kappa_group_count
                    write (biorthogonal_angular_unit) llistt(l), llistt(l)
                end if
            end do
            !
            !  Close the angular files
            !
            do l = 1, kappa_group_count
                close(kappa_coefficient_units(l), status="DELETE")
            end do
            !
            ncf0 = cumulative_csf_counts(iblk) + 1
        end do
        !
        !   Deallocate storage that is no longer required. This was
        !   allocated in lodcsl.
        !
        call dalloc (csf_subshell_occupations, "csf_subshell_occupations", "write_biorthogonal_mcp")
        call dalloc (csf_subshell_quantum_codes, "csf_subshell_quantum_codes", "write_biorthogonal_mcp")
        call dalloc (csf_coupling_quantum_codes, "csf_coupling_quantum_codes", "write_biorthogonal_mcp")
        return
    end subroutine write_biorthogonal_mcp
end module write_biorthogonal_mcp_interface
