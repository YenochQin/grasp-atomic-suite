!***********************************************************************
!                                                                      *
subroutine assemble_hamiltonian_block(symmetry_block_index, mpi_rank, mpi_process_count, total_contribution_count, &
    & total_secondary_contribution_count)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
    !      Spin-angular coefficiants have been put on memory               *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use grasp_parameters,    only: orbital_index_radix
    use array_memory_management, only: alloc, dalloc
    use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
        matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
        nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
    use mcp_file_metadata, only: max_mcp_file_multipole_rank
    use orbital_csf_basis, only: csf_count, orbital_count, orbital_kappa, orbital_two_j
    use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
        csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use rmcdhf_mcp_memory,     only: indx, coeff
    use csf_quantum_queries, only: iq
    use one_electron_radial_integral, only: rinti

    ! Procedure interfaces
    use direct_angular_coefficient_interface, only: direct_angular_coefficient

    ! Module dependencies
    use slater_radial_integral, only: slater

    ! Procedure interfaces
    use exchange_angular_coefficient_interface, only: exchange_angular_coefficient
    use load_mcp_records_into_memory_interface, only: load_mcp_records_into_memory
    use read_mcp_structure_into_memory_interface, only: read_mcp_structure_into_memory
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: symmetry_block_index
    integer, intent(in) :: mpi_rank, mpi_process_count
    integer(int64) :: total_contribution_count(20), total_secondary_contribution_count(20)

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix
    character(len=6), parameter :: routine_name = "SETHAM"

    ! Local declarations
    integer :: ncfpat, orbital_index, csf_index, itmp, idiag, km, k, k0, ib, nkjia, nkjib, kmin&
        , kmax, ndim, jblockt, ncft, ncoeff, lab, ncontr, iclmndum, i, &
        loc, nfile, id, ic
    integer(int64) :: lod      !  INTEGER*8 version of LOC
    real(real64) :: diaa, coef, f0aa, fkaa, f0ab, gkab, tegral
    logical :: set
    character :: mcp_channel_label*3
    !-----------------------------------------------
    !
    !     !*** needs EMT, IENDC.
    !     POINTER (PNTEMT,EMT(1))
    !     POINTER (PIENDC,IENDC(0:*))
    !     POINTER (PNIROW,IROW(1))
    !
    !     !*** needs ncf, nw
    !     POINTER (PNTRIQ,RIQDUMMY)
    !
    !
    !     POINTER (pncfpast, ncfpast(1))
    !     POINTER (pncminpast, ncminpast(1))
    !     POINTER (pnevecpast, nevecpast(1))
    !
    !     POINTER (PCOEFF,COEFF(1))
    !     !POINTER (PICLMN,ICLMNdum)
    !     POINTER (PINDX,INDX(1))
    !
    !-----------------------------------------------------------------------
    ncfpat = csf_offsets_by_block(symmetry_block_index)
    !=======================================================================
    !   Accumulate diagonal terms that do not require MCP coefficients
    !=======================================================================

    !=======================================================================
    !   Piece involving I(a,a) integrals
    !=======================================================================

    do orbital_index = 1, orbital_count
        set = .false.
        do csf_index = mpi_rank + 1, csf_count, mpi_process_count
            itmp = iq(orbital_index,csf_index + ncfpat)
            if (itmp <= 0) cycle
            !*** Occupation number not zero ...
            if (.not.set) then
                diaa = rinti(orbital_index,orbital_index,0)
                set = .true.
            end if
            ! IDIAG = IENDC(IR-1)+1
            idiag = stored_segment_end_offsets(csf_index)                    ! lower-triangle-by-rows mode
            hamiltonian_values(idiag) = hamiltonian_values(idiag) + itmp*diaa
        end do
    end do

    !=======================================================================
    !                    0
    !   Piece involving F (a,a) integrals
    !=======================================================================

    do orbital_index = 1, orbital_count
        set = .false.
        do csf_index = mpi_rank + 1, csf_count, mpi_process_count
            coef = direct_angular_coefficient(0,csf_index + ncfpat,orbital_index,orbital_index)
            if (coef == 0.0_real64) cycle
            !*** Angular coefficient not zero ...
            if (.not.set) then
                f0aa = slater(orbital_index,orbital_index,orbital_index,orbital_index,0)
                set = .true.
            end if
            ! IDIAG = IENDC(IR-1)+1
            idiag = stored_segment_end_offsets(csf_index)
            hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*f0aa
        end do
    end do

    !=======================================================================
    !                    k
    !   Piece involving F (a,a) integrals
    !=======================================================================

    km = 0
    k = 0
    do
        k = k + 2
        do orbital_index = 1, orbital_count
            k0 = orbital_two_j(orbital_index) - 1
            km = max0(k0,km)
            if (k > k0) cycle
            set = .false.

            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                coef = direct_angular_coefficient(k,csf_index + ncfpat,orbital_index,orbital_index)
                if (coef == 0.0_real64) cycle
                if (.not.set) then
                    fkaa = slater(orbital_index,orbital_index,orbital_index,orbital_index,k)
                    set = .true.
                end if
                ! IDIAG = IENDC(IR-1)+1
                idiag = stored_segment_end_offsets(csf_index)
                hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*fkaa
            end do
        end do
        if (k >= km) exit
    end do

    !=======================================================================
    !                    0
    !   Piece involving F (a,b) integrals
    !=======================================================================

    do orbital_index = 1, orbital_count - 1
        do ib = orbital_index + 1, orbital_count
            set = .false.
            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                coef = direct_angular_coefficient(0,csf_index + ncfpat,orbital_index,ib)
                if (coef == 0.0_real64) cycle
                if (.not.set) then
                    f0ab = slater(orbital_index,ib,orbital_index,ib,0)
                    set = .true.
                end if
                ! IDIAG = IENDC(IR-1)+1
                idiag = stored_segment_end_offsets(csf_index)
                hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*f0ab
            end do
        end do
    end do

    !=======================================================================
    !                    k
    !   Piece involving G (a,b) integrals
    !=======================================================================

    km = 0
    k = -1
    do
        k = k + 1
        do orbital_index = 1, orbital_count - 1
            nkjia = orbital_two_j(orbital_index)
            do ib = orbital_index + 1, orbital_count
                nkjib = orbital_two_j(ib)
                set = .false.
                if (orbital_kappa(orbital_index)*orbital_kappa(ib) > 0) then
                    kmin = abs((nkjia - nkjib)/2)
                else
                    kmin = abs((nkjia - nkjib)/2) + 1
                end if
                if (mod(k - kmin,2) /= 0) cycle

                kmax = (nkjia + nkjib)/2
                km = max0(kmax,km)
                if (k<kmin .or. k>kmax) cycle

                do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                    coef = exchange_angular_coefficient(k,csf_index + ncfpat,orbital_index,ib)
                    if (coef == 0.0_real64) cycle
                    if (.not.set) then
                        gkab = slater(orbital_index,ib,ib,orbital_index,k)
                        set = .true.
                    end if
                    ! IDIAG = IENDC(IR-1)+1
                    idiag = stored_segment_end_offsets(csf_index)
                    hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*gkab
                end do
            end do
        end do
        if (k >= km) exit
    end do

    !=======================================================================
    !   Local storage for reading mcpXXX files
    !=======================================================================

    ndim = 1
    call alloc (coeff, ndim, "COEFF", "SETHAM" )
    call alloc (indx, ndim, "INDX", "SETHAM")

    !=======================================================================
    !   Loop over non-zero labels which have non-zero elements
    !=======================================================================

    total_secondary_contribution_count(1) = total_secondary_contribution_count(1) + 1
    ! Cached labels are memory accesses; there is no IOSTAT result.
    call load_mcp_records_into_memory(31,total_secondary_contribution_count(1),lab, ncontr)
    do while(lab/=0 .or. ncontr/=0)

        !*** decode the label of I(ab)
        orbital_index = mod(lab,orbital_label_radix)
        ib = lab/orbital_label_radix

        !*** Compute radial integral I(ab)
        tegral = rinti(orbital_index,ib,0)

        ! Read column index, sparse matrix index, and coefficient
        ! for all contributions from this integral.
        if (ncontr > ndim) then
            call dalloc (coeff, "COEFF", "SETHAM")
            call dalloc (indx, "INDX", "SETHAM")
            ndim = ncontr
            call alloc (coeff, ndim, "COEFF", "SETHAM")
            call alloc (indx, ndim, "INDX", "SETHAM")
        end if
        call read_mcp_structure_into_memory(31,total_contribution_count(1),ncontr,iclmndum)
        total_contribution_count(1) = total_contribution_count(1) + ncontr

        !*** Store all the contributions from this integral
        do i = 1, ncontr
            loc = indx(i)
            lod = loc   !  (Convert type)
            if (lod > stored_matrix_element_count) then
                write (output_unit, *) "  Error in computing 1-e contribution"
                write (output_unit, *) "  LOC = ", loc , "  NELMNT = ", stored_matrix_element_count
                stop
            end if
            hamiltonian_values(loc) = hamiltonian_values(loc) + tegral*coeff(i)
        end do

        total_secondary_contribution_count(1) = total_secondary_contribution_count(1) + 1
        call load_mcp_records_into_memory(31,total_secondary_contribution_count(1),lab, ncontr)
    end do

    do nfile = 32, 32 + max_mcp_file_multipole_rank
        k = nfile - 32

        !=======================================================================
        !   Loop over non-zero labels which have non-zero elements
        !=======================================================================

        total_secondary_contribution_count(nfile-30) = total_secondary_contribution_count(nfile-30) + 1
        call load_mcp_records_into_memory(nfile,total_secondary_contribution_count(nfile-30),lab, ncontr)
        do while(lab/=0 .or. ncontr/=0)

            !                         k
            !*** decode the label of R (abcd)
            id = mod(lab,orbital_label_radix)
            lab = lab/orbital_label_radix
            ib = mod(lab,orbital_label_radix)
            lab = lab/orbital_label_radix
            ic = mod(lab,orbital_label_radix)
            orbital_index = lab/orbital_label_radix

            !*** Compute radial integral
            tegral = slater(orbital_index,ib,ic,id,k)

            ! Read column index, sparse matrix index, and coefficient
            ! for all contributions from this integral.
            if (ncontr > ndim) then
                call dalloc (coeff, "COEFF", "SETHAM")
                call dalloc (indx, "INDX", "SETHAM")
                ndim = ncontr
                call alloc (coeff, ndim, "COEFF", "SETHAM")
                call alloc (indx, ndim, "INDX", "SETHAM")
            end if
            call read_mcp_structure_into_memory(nfile,total_contribution_count(nfile-30),ncontr,iclmndum)
            total_contribution_count(nfile-30) = total_contribution_count(nfile-30) + ncontr

            !*** Store all the contributions from this integral
            do i = 1, ncontr
                loc = indx(i)
                lod = loc
                if (lod > stored_matrix_element_count) then   !! NELMENT (from Hmat_C) is INT*8
                    write (output_unit, *) "  Error in computing 2-e contribution"
                    write (output_unit, *) "  LOC= ", loc , "  NELMNT = ", stored_matrix_element_count
                    stop
                end if
                hamiltonian_values(loc) = hamiltonian_values(loc) + tegral*coeff(i)
            end do

            total_secondary_contribution_count(nfile-30) = total_secondary_contribution_count(nfile-30) + 1
            call load_mcp_records_into_memory(nfile,total_secondary_contribution_count(nfile-30),lab, ncontr)
        end do
    end do

    !=======================================================================
    !   Deallocate local storage
    !=======================================================================

    call dalloc (coeff, "COEFF", "SETHAM")
    call dalloc (indx, "INDX", "SETHAM")

    return
end subroutine assemble_hamiltonian_block
