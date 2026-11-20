!***********************************************************************
!                                                                      *
subroutine calculate_hyperfine_matrix_elements
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    use calculation_file_units, only: hyperfine_output_unit, summary_output_unit
    !                                                                      *
    !   This routine controls the main sequence of routine calls for the   *
    !   calculation  of the  hyperfine structure parameters.               *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, DRACAH, ISPAR, ITJPO,  *
    !                        TNSRJJ.                                       *
    !               [HFS92]: MATELT, RINT, RINTHF.                         *
    !                                                                      *
    !   Written by Per Jonsson and Farid A. Parpia                         *
    !                                                                      *
    !   Modified by Per Jonsson to evaluate g_j factors                    *
    !                                                                      *
    !                                         Last revision: 22 Oct 1999   *
    !   Modified by Gediminas Gaigalas for new spin-angular integration.   *
    !                                         Last revision:    Nov 2017   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals
    use array_memory_management, only: alloc, dalloc
    use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
    use atomic_calculation_state,           only: electron_mass_u, vacuum_speed_of_light_au, proton_mass_u, bohr_radius_cm, &
        & hartree_to_wavenumber_cm,     &
        speed_of_light_cm_per_s, b1=>electron_mass_u
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use first_order_csf_cutoffs, only: iccut
    use angular_symmetry_labels, only: right_justified_j_labels => jlbr, parity_labels => jlbp
    use isotope_properties,        only: sqn, dmomnm, qmomb,                  &
        hfsi=>sqn, hfsd=>dmomnm, hfsq=>qmomb
    use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0

    ! Procedure interfaces
    use hyperfine_radial_integral_interface, only: hyperfine_radial_integral

    ! Module dependencies
    use weighted_radial_overlap, only: rint

    ! Procedure interfaces
    use hyperfine_one_electron_element_interface, only: hyperfine_one_electron_element

    ! Module dependencies
    use integer_format, only: convrt
    use csf_quantum_queries, only: ispar, itjpo
    use one_electron_tensor_coefficients, only: oneparticlejj
    use racah_coefficients_scott_fischer, only: gracah1
    implicit none (type, external)

    ! Local constants
    real(real64), parameter :: cutoff = 1.0e-10_real64

    ! Local declarations
    integer :: ffmin, ffmax, ff, i, orbital_index, tensor_rank, operator_parity, csf_index, text_length, ir, itjpoc, itjpor&
        , idiff, ia, ket_orbital_index, k, kk, loc1, loc2, ii, jj, jjii, jb, ja, jjb, jja, iflag
    real(real64) :: shell_coefficients(max_orbitals)
    real(real64) :: rintme(2, max_orbitals, max_orbitals), amelt(2, max_orbitals, max_orbitals)
    real(real64) :: rintgj(max_orbitals, max_orbitals), rintdgj(max_orbitals, max_orbitals), &
        gjmelt(max_orbitals, max_orbitals), dgjmelt(max_orbitals, max_orbitals)
    !     .. Local pointer arrays
    real(real64), pointer :: hfc(:, :)
    real(real64), pointer :: gjc(:), dgjc(:)
    real(real64) ::  apart, gjpart, dgjpart, elemnt,  elemntgj, elemntdgj,&
        contr, contrgj, contrdgj, aumhz, barnau, dnmau, gfac, hfac, fj, &
        gja1, afa1, afa2, bfa1, bfa2, bfa3, gj, dgj, tilde1, &
        tilde2, factor1, factor2, rac1, rac2, hfselt1, hfselt2, evecic1, evecir1,&
        evecic2, evecir2
    character :: number_text*11
    !-----------------------------------------------
    !
    !
    !   Matrix elements smaller than CUTOFF are not accumulated
    !
    !
    !   Allocate storage for local arrays
    !
    call alloc (hfc, 5, nvec*nvec, "HFC", "HFS")
    call alloc (gjc, nvec*nvec, "GJC", "HFS")
    call alloc (dgjc, nvec*nvec, "DGLC", "HFS")
    !
    !   Initialise
    !
    hfc(:,:nvec*nvec) = 0.0_real64
    !
    gjc(:nvec*nvec) = 0.0_real64
    dgjc(:nvec*nvec) = 0.0_real64

    !
    !   Calculate and save the radial integrals and angular
    !   matrix elements for the two multipolarities
    !
    do tensor_rank = 1, 2
        do i = 1, orbital_count
            do orbital_index = 1, orbital_count
                if (tensor_rank == 1) then
                    rintme(tensor_rank,i,orbital_index) = hyperfine_radial_integral(i,orbital_index,-2)
                    rintgj(i,orbital_index) = hyperfine_radial_integral(i,orbital_index,1)
                else
                    rintme(tensor_rank,i,orbital_index) = rint(i,orbital_index,-3)
                    rintdgj(i,orbital_index) = rint(i,orbital_index,0)
                end if
                call hyperfine_one_electron_element (i, tensor_rank, orbital_index, apart, gjpart, dgjpart)
                amelt(tensor_rank,i,orbital_index) = apart
                if (tensor_rank /= 1) cycle
                gjmelt(i,orbital_index) = gjpart
                dgjmelt(i,orbital_index) = dgjpart
            end do
        end do
    end do
    !
    !   Set the parity of the one-body operators
    !
    operator_parity = 1
    !
    !   Sweep through the Hamiltonian matrix to determine the
    !   diagonal and off-diagonal hyperfine constants
    !
    do csf_index = 1, csf_count
        itjpoc = itjpo(csf_index)
        !
        !   Output IC on the screen to show how far the calculation has preceede
        !
        if (mod(csf_index,100) == 0) then
            call convrt (csf_index, number_text, text_length)
            write (output_unit, *) "Column "//number_text(1:text_length)//" complete;"
        end if
        !
        do ir = 1, csf_count
            !
            !   If LFORDR is .TRUE., a `first order' calculation is indicated;
            !   only the CSFs with serial numbers exceeding IC are treated specially
            !   only diagonal elements are evaluated for the `first order' CSFs
            !
            if (lfordr .and. csf_index>iccut .and. csf_index/=ir) cycle
            !
            itjpor = itjpo(ir)
            idiff = itjpoc - itjpor
            !
            !   Loop over the multipolarities
            !
            do tensor_rank = 1, 2
                !
                !   Initialise the accumulator
                !
                elemnt = 0.0_real64
                elemntgj = 0.0_real64
                elemntdgj = 0.0_real64
                !
                !   Consider 3 cases
                !               (k)
                !   (1) < J || T   || J >    , k = 1,2  and IR >= IC
                !
                !               (k)
                !   (2) < J || T   || J-1 >  , k = 1,2
                !
                !               (k)
                !   (3) < J || T   || J-2 >  , k = 2
                !
                if (.not.(idiff==0 .and. ir>=csf_index .or. idiff==2 .or. idiff==4&
                    .and. tensor_rank==2)) cycle
                !
                call oneparticlejj(tensor_rank,operator_parity,csf_index,ir,ia,ket_orbital_index,shell_coefficients)
                !GG               CALL TNSRJJ (KT, IPT, IC, IR, IA, IB, TSHELL)
                !
                !   Accumulate the contribution from the one-body operators;
                !
                if (ia /= 0) then
                    if (ia == ket_orbital_index) then
                        if (tensor_rank/=1 .or. idiff/=0) then
                            do ia = 1, orbital_count
                                if (abs(shell_coefficients(ia)) <= cutoff) cycle
                                elemnt = elemnt + amelt(tensor_rank,ia,ia)*rintme(tensor_rank,ia,ia)*&
                                    shell_coefficients(ia)
                                cycle
                            end do
                        else
                            do ia = 1, orbital_count
                                if (abs(shell_coefficients(ia)) <= cutoff) cycle
                                elemnt = elemnt + amelt(tensor_rank,ia,ia)*rintme(tensor_rank,ia,ia)*&
                                    shell_coefficients(ia)
                                elemntgj = elemntgj + gjmelt(ia,ia)*rintgj(ia,ia)*&
                                    shell_coefficients(ia)
                                elemntdgj = elemntdgj + dgjmelt(ia,ia)*rintdgj(ia,ia&
                                    )*shell_coefficients(ia)
                            end do
                        end if
                    else
                        if (abs(shell_coefficients(1)) > cutoff) then
                            elemnt = elemnt + amelt(tensor_rank,ia,ket_orbital_index)*rintme(tensor_rank,ia,ket_orbital_index)*&
                                shell_coefficients(1)
                            if (tensor_rank==1 .and. idiff==0) then
                                elemntgj = elemntgj + gjmelt(ia,ket_orbital_index)*rintgj(ia,ket_orbital_index)*&
                                    shell_coefficients(1)
                                elemntdgj = elemntdgj + dgjmelt(ia,ket_orbital_index)*rintdgj(ia,ket_orbital_index&
                                    )*shell_coefficients(1)
                            end if
                        end if
                    end if
                end if
                !
                !   Multiply with the configuration expansion coefficients and add the
                !   contributions from the matrix elements to obtain total contributions
                !
                if (abs(elemnt) <= cutoff .and. abs(elemntgj) <= cutoff .and. &
                    abs(elemntdgj) <= cutoff) cycle
                do k = 1, nvec
                    loc1 = (k - 1)*csf_count
                    evecic1 = csf_mixing_coefficients(csf_index + loc1)
                    evecir1 = csf_mixing_coefficients(ir + loc1)
                    do kk = 1, nvec
                        loc2 = (kk - 1)*csf_count
                        evecic2 = csf_mixing_coefficients(csf_index + loc2)
                        evecir2 = csf_mixing_coefficients(ir + loc2)
                        if (idiff==0 .and. ir/=csf_index) then
                            contr = elemnt*(evecic1*evecir2 + evecir1*evecic2)
                            contrgj = elemntgj*(evecic1*evecir2 + evecir1*evecic2)
                            contrdgj = elemntdgj*(evecic1*evecir2 + evecir1*evecic2)
                        else
                            contr = elemnt*evecic1*evecir2
                            contrgj = elemntgj*evecic1*evecir2
                            contrdgj = elemntdgj*evecic1*evecir2
                        end if
                        !
                        !   Magnetic dipole and the two operators of the g_j factor
                        !
                        if (tensor_rank == 1) then
                            if (idiff == 0) then
                                hfc(1,nvec*(k-1)+kk) = hfc(1,nvec*(k-1)+kk) + contr
                                gjc(nvec*(k-1)+kk) = gjc(nvec*(k-1)+kk) + contrgj
                                dgjc(nvec*(k-1)+kk) = dgjc(nvec*(k-1)+kk) + contrdgj
                            else if (itjpoc - itjpor == 2) then
                                hfc(2,nvec*(k-1)+kk) = hfc(2,nvec*(k-1)+kk) + contr
                            end if
                        !
                        !   Electric quadrupole
                        !
                        else if (tensor_rank == 2) then
                            select case (idiff)
                            case (0)
                                hfc(3,nvec*(k-1)+kk) = hfc(3,nvec*(k-1)+kk) + contr
                            case (2)
                                hfc(4,nvec*(k-1)+kk) = hfc(4,nvec*(k-1)+kk) + contr
                            case (4)
                                hfc(5,nvec*(k-1)+kk) = hfc(5,nvec*(k-1)+kk) + contr
                            end select
                        end if
                    end do
                end do
                !
            end do
            !
        end do
    end do
    !
    !   These are the conversion factors to obtain the hyperfine
    !   constants in MHz
    !
    aumhz = hartree_to_wavenumber_cm*speed_of_light_cm_per_s*1.0e-6_real64
    barnau = 1.0e-24_real64/bohr_radius_cm**2
    dnmau = b1/(2.0_real64*vacuum_speed_of_light_au*proton_mass_u)
    !
    gfac = aumhz*dnmau*hfsd/hfsi
    hfac = aumhz*2.0_real64*hfsq*barnau
    !
    !   Output the hyperfine interaction constants
    !
    write (summary_output_unit, '(/,/," Interaction constants:"/,/," Level1  J Parity  Level2  J Parity",8x,"A (MHz)",13x,&
        &"B (MHz)"/)')
    write (hyperfine_output_unit, '(/,/," Interaction constants:"/,/," Level1  J Parity ",8x,"A (MHz)",13x,"B (MHz)",13x,&
        &"g_J",14x,&
        &"delta g_J",11x,"total g_J"/)')
    !
    do i = 1, nvec
        do ii = 1, nvec
            !
            jj = iatjpo(i)
            jjii = iatjpo(ii)
            !
            if (.not.(jj==jjii .and. jjii>1 .or. jj>jjii)) cycle
            !
            fj = 0.5_real64*dble(jj - 1)
            !
            gja1 = sqrt(1.0_real64/(fj*(fj + 1.0_real64)))

            afa1 = gfac*gja1

            if (jj == 2) then
                afa2 = 0.0_real64
            else
                afa2 = gfac*sqrt(1.0_real64/(fj*(2.0_real64*fj - 1.0_real64)))
            end if
            bfa1 = hfac*sqrt((fj*(2.0_real64*fj - 1.0_real64))/((fj + 1.0_real64)*(2.0_real64*&
                fj + 3.0_real64)))
            if (jj == 2) then
                bfa2 = 0.0_real64
            else
                bfa2 = 0.25_real64*hfac*sqrt((fj*(fj - 1.0_real64))/((fj + 1.0_real64)*(&
                    2.0_real64*fj - 1.0_real64)))
            end if
            if (jj == 4) then
                bfa3 = 0.0_real64
            else
                bfa3 = 0.125_real64*hfac*sqrt((fj*(fj - 1.0_real64)*(2.0_real64*fj - 1.0_real64&
                    ))/(2.0_real64*fj - 3.0_real64))
            end if
            !
            !   Diagonal (J,J) A and B factors
            !
            select case (jj - jjii)
            case (0)
                if (i <= ii) then
                    !                GJ = CVAC*GJA1*GJC(NVEC*(I-1)+II)
                    !                DGJ = 0.001160D0*GJA1*DGJC(NVEC*(I-1)+II)

                    write (summary_output_unit, '(1x,1i3,5x,2a4,2x,1i3,5x,2a4,1p,2d20.10)') ivec(i), right_justified_j_labels(jj), &
                        parity_labels((iaspar(i)+3)/2), &
                        ivec(ii), right_justified_j_labels(jjii), parity_labels((iaspar(ii)+3)/2), afa1*hfc(1,&
                        nvec*(i-1)+ii), bfa1*hfc(3,nvec*(i-1)+ii)
                    !
                    !   Output diagonal hfs and g_j factors to file <name>.h or <name>.ch
                    !
                    if (i == ii) then
                        gj = vacuum_speed_of_light_au*gja1*gjc(nvec*(i-1)+ii)
                        dgj = 0.001160_real64*gja1*dgjc(nvec*(i-1)+ii)
                        write (hyperfine_output_unit, '(1x,1i3,5x,2a4,1p,5d20.10)') ivec(i), right_justified_j_labels(jj), &
                            & parity_labels((iaspar(i)+3)/2), &
                            afa1*hfc(1,nvec*(i-1)+ii), bfa1*hfc(3,nvec*(i-1)+ii), &
                            gj, dgj, gj + dgj
                    end if
                end if
            !
            !   Off diagonal (J,J-1) A and B factors
            !
            case (2)
                write (summary_output_unit, '(1x,1i3,5x,2a4,2x,1i3,5x,2a4,1p,2d20.10)') &
                    ivec(i), right_justified_j_labels(jj), parity_labels((iaspar(i)+3)/2), ivec(&
                    ii), right_justified_j_labels(jjii), parity_labels((iaspar(ii)+3)/2), afa2*hfc(2,nvec*(i-1&
                    )+ii), bfa2*hfc(4,nvec*(i-1)+ii)
            !
            !   Off diagonal (J,J-2) B factor
            !
            case (4)
                write (summary_output_unit, '(1x,1i3,5x,2a4,2x,1i3,5x,2a4,1p,2d20.10)') &
                    ivec(i), right_justified_j_labels(jj), parity_labels((iaspar(i)+3)/2), ivec(&
                    ii), right_justified_j_labels(jjii), parity_labels((iaspar(ii)+3)/2), 0.0_real64, bfa3*hfc(5,&
                    nvec*(i-1)+ii)
            end select
            !
        end do
    end do
    !
    !   These are the factors needed to obtain the F-dependent hyperfine
    !   matrix elements in Hartrees
    !
    tilde1 = sqrt((hfsi + 1.0_real64)/hfsi)*hfsd*0.5_real64*b1/(proton_mass_u*vacuum_speed_of_light_au)
    if (hfsi > 0.6_real64) then
        tilde2 = sqrt((hfsi + 1.0_real64)*(2.0_real64*hfsi + 3.0_real64)/((2.0_real64*hfsi - &
            1.0_real64)*hfsi))
        tilde2 = tilde2*hfsq*0.5_real64*barnau
    else
        tilde2 = 0.0_real64
    end if
    !
    !     II = 2*I
    !
    ii = nint(2.0_real64*hfsi)
    !
    !   Calculate the F-dependent matrix elemnts
    !
    !   Loop over the states
    !
    do jb = 1, nvec
        do ja = 1, nvec
            jjb = iatjpo(jb) - 1
            jja = iatjpo(ja) - 1
            if (.not.(jja==jjb .and. jja>0 .or. jjb>jja)) cycle
            !
            !   Determine the possible F quantum numbers for the matrix element
            !
            ffmin = max(abs(jja - ii),abs(jjb - ii))
            ffmax = min(jja + ii,jjb + ii)
            !
            !   Loop over the possible F quantum numbers
            !
            iflag = 0
            do ff = ffmin, ffmax, 2
                !
                !   Phase factor
                !
                if (mod((ii + jja - ff)/2,2) == 1) then
                    factor1 = -1.0_real64
                else
                    factor1 = 1.0_real64
                end if
                factor2 = factor1*sqrt((dble(jjb) + 1.0_real64)*(dble(ii) + 1.0_real64))
                !
                !   Determine the Racah W coefficients.
                !
                !GG               CALL DRACAH (II, JJA, II, JJB, FF, 2, RAC1)
                call gracah1 (ii,jja,ii,jjb,ff,2,rac1)
                !GG               CALL DRACAH (II, JJA, II, JJB, FF, 4, RAC2)
                call gracah1 (ii,jja,ii,jjb,ff,4,rac2)
                !
                !   Obtain and output matrix elements for J,J
                !
                if (jja - jjb == 0) then
                    hfselt1 = factor2*rac1*hfc(1,nvec*(jb-1)+ja)*tilde1
                    hfselt2 = factor2*rac2*hfc(3,nvec*(jb-1)+ja)*tilde2
                    if (abs(hfselt1 + hfselt2) > cutoff*10.0e-5_real64) then
                        if (iflag == 0) write (summary_output_unit, '(/,/," Matrix elements:"/,/,&
                            &" Level1  J Parity  Level2  J Parity    F ",4x,"Matrix element (a.u.)"/)')
                        iflag = 1
                        write (summary_output_unit, '(1x,1i3,5x,2a4,2x,1i3,5x,2a4,1x,a4,4x,1p,1d20.10)') &
                            ivec(jb), right_justified_j_labels(jjb+1), parity_labels((iaspar(jb)+3)&
                            /2), ivec(ja), right_justified_j_labels(jja+1), parity_labels((iaspar(ja)+3)/2), &
                            right_justified_j_labels(ff+1), hfselt1 + hfselt2
                    end if
                !
                !   Obtain and output matrix elements for J,J-1
                !
                else if (abs(jja - jjb) == 2) then
                    hfselt1 = factor2*rac1*hfc(2,nvec*(jb-1)+ja)*tilde1
                    hfselt2 = factor2*rac2*hfc(4,nvec*(jb-1)+ja)*tilde2
                    if (abs(hfselt1 + hfselt2) > cutoff*10.0e-5_real64) then
                        if (iflag == 0) write (summary_output_unit, '(/,/," Matrix elements:"/,/,&
                            &" Level1  J Parity  Level2  J Parity    F ",4x,"Matrix element (a.u.)"/)')
                        iflag = 1
                        write (summary_output_unit, '(1x,1i3,5x,2a4,2x,1i3,5x,2a4,1x,a4,4x,1p,1d20.10)') &
                            ivec(jb), right_justified_j_labels(jjb+1), parity_labels((iaspar(jb)+3)&
                            /2), ivec(ja), right_justified_j_labels(jja+1), parity_labels((iaspar(ja)+3)/2), &
                            right_justified_j_labels(ff+1), hfselt1 + hfselt2
                    end if
                !
                !   Obtain and output matrix elements for J,J-2
                !
                else if (abs(jja - jjb) == 4) then
                    hfselt1 = 0.0_real64
                    hfselt2 = factor2*rac2*hfc(5,nvec*(jb-1)+ja)*tilde2
                    if (abs(hfselt1 + hfselt2) > cutoff*10.0e-5_real64) then
                        if (iflag == 0) write (summary_output_unit, '(/,/," Matrix elements:"/,/,&
                            &" Level1  J Parity  Level2  J Parity    F ",4x,"Matrix element (a.u.)"/)')
                        iflag = 1
                        write (summary_output_unit, '(1x,1i3,5x,2a4,2x,1i3,5x,2a4,1x,a4,4x,1p,1d20.10)') &
                            ivec(jb), right_justified_j_labels(jjb+1), parity_labels((iaspar(jb)+3)&
                            /2), ivec(ja), right_justified_j_labels(jja+1), parity_labels((iaspar(ja)+3)/2), &
                            right_justified_j_labels(ff+1), hfselt1 + hfselt2
                    end if
                else
                    hfselt1 = 0.0_real64
                    hfselt2 = 0.0_real64
                end if
            end do
        end do
    end do
    !
    call dalloc (hfc, "HFC", "HFS")
    return
    !
    return
    !
end subroutine calculate_hyperfine_matrix_elements
