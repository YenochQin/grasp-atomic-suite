!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
subroutine compute_mcp_angular_coefficients_mpi (nb, is_restart, output_file_prefix)
    ! Module dependencies
    use mcp_generation_options, only: mcp_temporary_units
    use calculation_file_units, only: summary_output_unit
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This routine controls the computation  and storage of the values   *
    !   and all indices of the angular coefficients                        *
    !                                                                      *
    !                                       k                              *
    !                   T  (ab)            V  (abcd)                       *
    !                    rs                 rs                             *
    !                                                                      *
    !   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
    !   c and d are orbital sequence numbers.  r and s are configuration   *
    !   state function indices.                                            *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, RKCO_GG,       *
    !                        TNSRJJ.                                       *
    !               [GENMCP]: FNDBEG, SETSDA, SORT.                        *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
    !   Modified by C. Froese Fischer for block computation.               *
    !   Modified by G. Gaigalas and J. Bieron for new spin-angular         *
    !   integration                                        01 April 2012   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use iso_fortran_env, only: int64
    use rangular_profile, only: stage_clock, record_stage
    use grasp_parameters,   only:  max_orbitals, orbital_index_radix
    use array_memory_management, only: alloc, dalloc
    use mpi_state,           only:  mpi_rank, mpi_process_count

    ! Procedure interfaces
    use write_mcp_sparsity_summary_mpi_interface, only: write_mcp_sparsity_summary_mpi
    use locate_mcp_restart_position_interface, only: locate_mcp_restart_position

    ! Module dependencies
    use angular_coefficient_buffer, only: alcbuf
    use two_electron_angular_coefficients, only: rkco_gg

    ! Procedure interfaces
    use write_mcp_sparsity_interface, only: write_mcp_sparsity
    use sort_mcp_coefficients_interface, only: sort_mcp_coefficients

    ! Module dependencies
    use angular_coefficient_buffer_state,  only: angular_coefficient_count, angular_coefficient_labels, angular_coefficients
    use debug_output_settings,   only: ldbpa
    use calculation_input_defaults, only: ndef
    use csf_block_cutoffs,    only: iccut
    use mcp_generation_options,     only: max_coulomb_multipole_rank, diagonal_only, use_first_order_restriction
    use orbital_csf_basis,     only: csf_count, csf_subshell_occupations
    use csf_coupling_data,    only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
    use one_electron_scalar_coefficients, only: onescalar
    use closed_shell_mcp, only: cord
    implicit none (type, external)

    ! Arguments
    integer  :: nb
    !     INTEGER  :: MYID
    !     INTEGER  :: NPROCS
    logical , intent(in) :: is_restart
    character(len=*), intent(in) :: output_file_prefix

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix
    integer, parameter :: keysq = orbital_index_radix*orbital_index_radix
    !GG      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
    !GG      REAL(DOUBLE), PARAMETER :: CUTOFF = 0.0D0
    !cjb  REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-08
    real(real64), parameter :: cutoff = 1.0e-10_real64

    ! Local declarations
    integer, pointer :: llistv(:)
    integer :: first_bra_csf_index,first_ket_csf_index,npos,llistt,k,csf_index,ket_csf_index,bra_orbital_index, &
        & ket_orbital_index,lab,i,ic, &
        id, nswap, iswap, lac, lbd, ntgi
    real(real64) :: shell_coefficients(max_orbitals)
    real(real64) :: vcoeff
    logical :: f0int, lincr
    integer(int64) :: stage_start_count
    !-----------------------------------------------
    !
    !cjb           PRINT *, ' mcpmpi_gg.f90 FHEAD '
    !cjb           PRINT *, ' mcpmpi_gg.f90 FHEAD = ', FHEAD
    !
    !   Allocate storage to the array that stores the lengths
    !   of the lists of V coefficients of each multipolarity
    !
    stage_start_count = stage_clock()
    call alloc(llistv, 0, max_coulomb_multipole_rank, "LLISTV", "MCPMPI")
    !
    !   If this is a restart, determine the pair of CSFs with which
    !   the computation should begin.
    !   If not restart, initialize the counters for the total number of
    !   T coefficients and V coefficients of each multipolarity.
    !
    if (is_restart) then
        print *, "Not ready for RESTART"
        stop
        call locate_mcp_restart_position (first_bra_csf_index, first_ket_csf_index, npos, llistt, llistv)
    else
        first_bra_csf_index = 1 + mpi_rank
        npos = 0
        llistt = 0
        do k = 0, max_coulomb_multipole_rank
            llistv(k) = 0
        end do
    end if
    !
    !   Allocate storage for the arrays in BUFFER
    !
    call alcbuf (1)
    !
    !   Write header to  .dbg  file if appropriate
    !
    if (ldbpa(2) .or. ldbpa(3)) then
        call ensure_debug_output()
        write (angular_debug_unit,'(/"From MCPMPI:")')
    end if
    !
    !   JA and JB respectively refer to the initial and final states
    !   in the list of NCF configurations
    !
    call record_stage("coefficient_buffer_allocation",nb,stage_start_count)
    stage_start_count = stage_clock()
    do csf_index = first_bra_csf_index, csf_count, mpi_process_count
        first_ket_csf_index = 1
        do ket_csf_index = first_ket_csf_index, csf_index
            if (diagonal_only .or. (use_first_order_restriction .and. (ket_csf_index > iccut(nb))) ) then
                if (ket_csf_index /= csf_index) cycle
            end if
            !
            !   LINCR is .TRUE. if npos is to be incremented by 1; there
            !   is always a diagonal element in each column
            !
            if (ket_csf_index /= csf_index) then
                lincr = .true.
            else
                npos = npos + 1
                lincr = .false.
            end if
            if (ket_csf_index /= csf_index) then
                !   Compute T coefficients
                call onescalar(csf_index,ket_csf_index,bra_orbital_index,ket_orbital_index,shell_coefficients)
                if (bra_orbital_index /= 0 .and. bra_orbital_index /= ket_orbital_index .and.       &
                    abs (shell_coefficients(1)) > cutoff) then
                    npos = npos + 1
                    lincr = .false.
                    llistt = llistt + 1
                    lab = min (bra_orbital_index,ket_orbital_index) * orbital_label_radix + max (bra_orbital_index, &
                        & ket_orbital_index)
                    write (mcp_temporary_units(31)) csf_index, npos, lab, shell_coefficients(1)
                end if
            end if
            !
            !   Call the MCP package to generate V coefficients; ac and bd
            !   are the density pairs. NVCOEF is initialized here but changed
            !   in /rkco/cor[d] via COMMON.
            !
            angular_coefficient_count = 0
            !            CALL RKCO_GG (JA, JB, 0, 1)
            call rkco_gg (csf_index, ket_csf_index, cord, 0, 1)
            do i = 1, angular_coefficient_count
                vcoeff = angular_coefficients(i)
                if (abs (vcoeff) > cutoff) then
                    bra_orbital_index = angular_coefficient_labels(1,i)
                    ket_orbital_index = angular_coefficient_labels(2,i)
                    ic = angular_coefficient_labels(3,i)
                    id = angular_coefficient_labels(4,i)
                    k  = angular_coefficient_labels(5,i)
                    f0int=(k==0).and.(bra_orbital_index==ic).and.(ket_orbital_index==id)
                    if (.not. f0int) then
                        ! Swap index to make sure IA <= IC, IB <= ID and record the number
                        ! of swaps
                        nswap = 0
                        if (bra_orbital_index > ic) then
                            iswap = ic
                            ic = bra_orbital_index
                            bra_orbital_index = iswap
                            nswap = nswap + 1
                        end if
                        if (ket_orbital_index > id) then
                            iswap = id
                            id = ket_orbital_index
                            ket_orbital_index = iswap
                            nswap = nswap + 1
                        end if
                        if (lincr) then
                            npos = npos + 1
                            lincr = .false.
                        end if
                        llistv(k) = llistv(k) + 1
                        lac = bra_orbital_index * orbital_label_radix + ic
                        lbd = ket_orbital_index * orbital_label_radix + id
                        if (lac < lbd) then
                            lab = lac * keysq + lbd
                        else
                            lab = lbd * keysq + lac
                        end if
                        write (mcp_temporary_units(32+k)) csf_index, npos, lab, vcoeff
                    end if
                end if
            end do
            !
            !   All angular coefficients for this pair of CSFs have been
            !   generated; update file 30
            !
            if ((ket_csf_index == csf_index) .or. (.not. lincr))   &
                write (mcp_temporary_units(30)) csf_index, ket_csf_index, npos
        end do
        if (mod(csf_index-1,100) == 0) print "(A4,I7,A9,I10,A9,I3,A8,I3)",  &
            "Row ",csf_index," nnonz = ",npos," block = ",nb," myid = ",mpi_rank
    end do
    call record_stage("pair_generation_and_buffering",nb,stage_start_count, &
        int(llistt,int64)+sum(int(llistv,int64)),storage="disk")
    stage_start_count = stage_clock()
    !
    !   Deallocate storage that is no longer required
    !
    call dalloc (csf_subshell_occupations, "IQA", "MCPMPI")
    call dalloc (csf_subshell_quantum_codes, "JQSA", "MCPMPI")
    call dalloc (csf_coupling_quantum_codes, "JCUPA", "MCPMPI")
    call alcbuf (3)
    call record_stage("csf_and_pair_buffer_release",nb,stage_start_count)
    !
    !   Write out a report for this run
    !
    if (ndef /= 0 .and. mpi_rank == 0) then
        write (summary_output_unit,*)
        write (summary_output_unit,*) llistt, " T coefficients generated;"
        do k = 0, max_coulomb_multipole_rank
            write (summary_output_unit,*) llistv(k)," V(k=",k,") coefficients generated;"
        end do
    end if
    !
    !   Set up sparse structure definition arrays in file 30
    !

    stage_start_count = stage_clock()
    rewind (mcp_temporary_units(30))
    call write_mcp_sparsity (write_mcp_sparsity_summary_mpi, npos, ldbpa(4), nb, mpi_rank, mpi_process_count, output_file_prefix)
    close (mcp_temporary_units(30), status = "DELETE")
    call record_stage("sparse_structure_and_output",nb,stage_start_count,int(npos,int64),storage="disk")
    !
    !   Sort MCP coefficients into integral-based lists
    !
    !   1) the T coefficient
    stage_start_count = stage_clock()
    rewind (mcp_temporary_units(31))
    call sort_mcp_coefficients (31, (llistt), ntgi, ldbpa(2), nb, output_file_prefix)
    if (ndef /= 0 .and. mpi_rank == 0) then
        write (summary_output_unit,*) ntgi, " I(ab) integrals;"
    end if
    close (mcp_temporary_units(31), status = "DELETE")
    call record_stage("sort_and_output",nb,stage_start_count,int(llistt,int64),31,"disk")
    !   2) the V coefficient
    do k = 0, max_coulomb_multipole_rank
        stage_start_count = stage_clock()
        rewind (mcp_temporary_units(k+32))
        call sort_mcp_coefficients (k+32, (llistv(k)), ntgi, ldbpa(3), nb, output_file_prefix)
        if (ndef /= 0 .and. mpi_rank == 0) then
            write (summary_output_unit,*) " k = ", k, ": ", ntgi, " Slater integrals;"
        end if
        close (mcp_temporary_units(k+32), status = "DELETE")
        call record_stage("sort_and_output",nb,stage_start_count,int(llistv(k),int64),k+32,"disk")
    end do
    call dalloc(llistv, "LLISTV", "MCPMPI")
    return
end subroutine compute_mcp_angular_coefficients_mpi
