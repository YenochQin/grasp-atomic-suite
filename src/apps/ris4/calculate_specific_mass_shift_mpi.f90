!***********************************************************************
!                                                                      *
subroutine calculate_specific_mass_shift_mpi(doit,vint,vint2)
    ! Module dependencies
    use calculation_file_units, only: two_body_cache_unit
    !                                                                      *
    !   MPI version of SMSNEW. Each rank handles a strided subset of the   *
    !   outer IC loop and can optionally write a per-rank angular-data      *
    !   fragment that the root process later merges into name.ITB.          *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: orbital_index_radix, max_orbitals
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count
    use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
        angular_coefficient_labels, angular_coefficients
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use isotope_shift_options, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
    use mpi_state,            only: mpi_rank, mpi_process_count
    use angular_coefficient_buffer, only: alcbuf
    use csf_quantum_queries, only: ispar, itjpo
    use integer_format, only: convrt
    use two_electron_angular_coefficients, only: rkco_gg
    use closed_shell_mcp, only: cord
    implicit none (type, external)

    ! Arguments
    real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
    integer, intent(in) :: doit

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix

    ! Local declarations
    character(len=11) :: number_text
    real(real64) :: evpair(nvec)
    real(real64) :: contri, contrik1, vcoeff, radial_k1, radial_sms,  &
        pair_scale, term_k1, term_sms
    integer :: j,k,loc,csf_index,ir,itjpoc,include_closed_shell_terms,text_length, iia, iib, iic, iid
    integer :: orbital_index, occ_ic, occ_ir, ndiff

    include_closed_shell_terms = 1
    call alcbuf (1)

    do csf_index = mpi_rank + 1, csf_count, mpi_process_count
        if (mpi_rank == 0 .and. mod(csf_index,100) == 1) then
            call convrt (csf_index,number_text,text_length)
            print *, "Column "//number_text(1:text_length)//" complete on root;"
        end if

        itjpoc = itjpo (csf_index)
        do ir = csf_index,csf_count
            if (itjpo(ir) == itjpoc) then
                if (ispar(csf_index) /= ispar(ir)) cycle
                ndiff = 0
                do orbital_index = 1,orbital_count
                    occ_ic = csf_subshell_occupations(orbital_index,csf_index)
                    occ_ir = csf_subshell_occupations(orbital_index,ir)
                    if (occ_ic == occ_ir) cycle
                    if (iabs(occ_ic - occ_ir) > 2) then
                        ndiff = 5
                        exit
                    end if
                    ndiff = ndiff + iabs(occ_ic - occ_ir)
                    if (ndiff > 4) exit
                end do
                if (ndiff > 4 .or. mod(ndiff,2) /= 0) cycle

                angular_coefficient_count = 0
                call rkco_gg (csf_index, ir, cord, include_closed_shell_terms, 1)
                if (angular_coefficient_count == 0) cycle
                if (ir == csf_index) then
                    pair_scale = 1.0_real64
                else
                    pair_scale = 2.0_real64
                end if
                do j = 1,nvec
                    loc = (j-1)*csf_count
                    evpair(j) = pair_scale * csf_mixing_coefficients(csf_index+loc) * csf_mixing_coefficients(ir+loc)
                end do
                do k = 1,angular_coefficient_count
                    vcoeff = angular_coefficients(k)
                    if (abs (vcoeff) > cutoff) then
                        if (angular_coefficient_labels(5,k) == 1) then
                            iia = angular_coefficient_labels(1,k)
                            iib = angular_coefficient_labels(2,k)
                            iic = angular_coefficient_labels(3,k)
                            iid = angular_coefficient_labels(4,k)
                            radial_k1 = vint(iia,iic) * vint(iib,iid)
                            radial_sms = (vint2(iia,iic) * vint(iib,iid) +   &
                                vint2(iib,iid) * vint(iia,iic))/2.0_real64
                            if (radial_k1 == 0.0_real64 .and.                     &
                                radial_sms == 0.0_real64) cycle
                            term_k1 = -vcoeff * radial_k1
                            term_sms = -vcoeff * radial_sms
                            if(doit==1) write(two_body_cache_unit) csf_index,ir
                            if(doit==1) then
                                write(two_body_cache_unit) vcoeff, ((iia*orbital_label_radix + iic)*orbital_label_radix+iib)* &
                                    & orbital_label_radix+iid
                            end if
                            do j = 1,nvec
                                contrik1 = evpair(j) * term_k1
                                contri = evpair(j) * term_sms
                                smsc1(j) = smsc1(j) + contrik1
                                smsc2(j) = smsc2(j) + contri
                            end do
                        end if
                    end if
                end do
            end if
        end do
    end do

    call alcbuf (3)
    return
end subroutine calculate_specific_mass_shift_mpi
