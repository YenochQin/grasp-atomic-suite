!***********************************************************************
!                                                                      *
subroutine densnew(doit,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
!       SUBROUTINE DENSNEW(DOIT,DINT1,DINT2,DINT3,                      &
!                          DINT4,DINT5,DINT6,DINT7,                     &
!                          DINT1VEC,DENS1VEC)
!                                                                      *
!   IF angular coefficients must be calculated                         *
!   This routine controls combines the radial and angular parts for the*
!   calculation of the NMS parameter, the electron density at the      *
!   origin and radial expectation values.
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, GETYN                         *
!                        ITJPO, ONESCALAR                              *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!   Modified by C. Naz\'e  Feb. 2012                                   *
!   Modified by J. Ekman   Nov. 2013                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: orbital_index_radix, max_orbitals, max_radial_points
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use atomic_calculation_state, only: max_radial_points
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use first_order_csf_cutoffs, only: iccut
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use angular_symmetry_labels, only: left_justified_j_labels => jlbl, parity_labels => jlbp, right_justified_j_labels => jlbr
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count, orbital_kappa
      use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
          nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use teilst_c, only: ndtea, ntei, pindte, pvalte, firstt
      use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
          angular_coefficient_labels, angular_coefficients
      use ris_c, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use selected_level_indices,           only : nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use angular_coefficient_buffer, only: alcbuf
      use integer_format, only: convrt
      use yes_no_input, only: getyn
      use csf_quantum_queries, only: itjpo
      use one_electron_scalar_coefficients, only: onescalar
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), &
                                                        dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                                        dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals), &
                                                        dint7(max_orbitals, max_orbitals)
      integer, intent(in) :: doit
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: shell_coefficients(max_orbitals), tshell_s(max_orbitals)
      real(real64), pointer :: emt1(:), emt2(:), emt3(:), emt4(:),   &
                                             emt5(:), emt6(:)
      real(real64) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5,     &
                      elemnt6, elemnt7, contri1, contri2, contri3,     &
                      contri4, contri5, contri6, contri7
      logical :: vsh, nucde, smssh, yes
      character :: number_text*11, ck*2
      integer :: ia_s(max_orbitals)
      integer :: ka, iopar, include_closed_shell_terms, csf_index, text_length, itjpoc, ir, orbital_index, ket_orbital_index, i, j
      integer :: l, loc, ncontr, lab, ii
!-----------------------------------------------
!
! DOIT: IF DOIT=1 angular coefficients will be stored after creation
! DINT1 contain the density
! DINT2 contain the uncorrected NMS parameter: K^1_NMS
! DINT3 contain the expect. value <r>
! DINT4 contain the expect. value <r2>
! DINT5 contain the expect. value <r-1>
! DINT6 contain the expect. value <r-2>
! DINT7 contain the sum of NMS parameters: K^1_NMS+K^2_NMS+K^3_NMS
!
!
!   Set the rank (zero) and parity (even) for the one-particle
!   coefficients
!
      ka = 0
      iopar = 1
      include_closed_shell_terms = 1
!
!   Allocate storage for the arrays in BUFFER
!
      call alcbuf (1)
!
!   Sweep through the Hamiltonian matrix to determine the
!   sms parameter
!
     do 13 csf_index = 1,csf_count
!
!   Output IC on the screen to show how far the calculation has preceede
!
        call convrt (csf_index,number_text,text_length)
        if (mod(csf_index,100)==0) then
          print *, "Column "//number_text(1:text_length)//" complete;"
        end if
!
        itjpoc = itjpo (csf_index)

        do 12 ir = csf_index,csf_count
!
!   Matrix elements are diagonal in J
!
          if (itjpo(ir) == itjpoc) then
!
!   Initialise the accumulator
!
            elemnt1 = 0.0_real64
            elemnt2 = 0.0_real64
            elemnt3 = 0.0_real64
            elemnt4 = 0.0_real64
            elemnt5 = 0.0_real64
            elemnt6 = 0.0_real64
            elemnt7 = 0.0_real64
!
!   Call the MCT package to compute T coefficients
!
           call onescalar(csf_index,ir,orbital_index,ket_orbital_index,shell_coefficients)
!GG            CALL TNSRJJ (KA,IOPAR,IC,IR,IA,IB,TSHELL)
            if (orbital_index /= 0) then
              if (orbital_index == ket_orbital_index) then
                    ncontr = 0
                do 8 orbital_index = 1,orbital_count
                  if (abs (shell_coefficients(orbital_index)) > cutoff) then
                    ncontr = ncontr + 1
                    tshell_s(ncontr) = shell_coefficients(orbital_index)
                    ia_s(ncontr) = orbital_index
                    elemnt1 = elemnt1 + dint1(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                    elemnt2 = elemnt2 + dint2(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                    elemnt3 = elemnt3 + dint3(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                    elemnt4 = elemnt4 + dint4(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                    elemnt5 = elemnt5 + dint5(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                    elemnt6 = elemnt6 + dint6(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                    elemnt7 = elemnt7 + dint7(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                  end if
    8           continue
                if (doit==1) write(50) csf_index,ir,ncontr
                do i = 1,ncontr
                  lab = ia_s(i)*(orbital_label_radix + 1)
                  if (doit==1) write(50) tshell_s(i),lab
                end do
              else
                if (abs (shell_coefficients(1)) > cutoff) then
                  if (orbital_kappa(orbital_index)==orbital_kappa(ket_orbital_index)) then
                    if (doit==1) write(50) csf_index,ir,1
                    lab = orbital_index*orbital_label_radix + ket_orbital_index
                    if (doit==1) write(50) shell_coefficients(1),lab
                    elemnt1 = elemnt1 + dint1(orbital_index,ket_orbital_index)*shell_coefficients(1)
                    elemnt2 = elemnt2 + dint2(orbital_index,ket_orbital_index)*shell_coefficients(1)
                    elemnt3 = elemnt3 + dint3(orbital_index,ket_orbital_index)*shell_coefficients(1)
                    elemnt4 = elemnt4 + dint4(orbital_index,ket_orbital_index)*shell_coefficients(1)
                    elemnt5 = elemnt5 + dint5(orbital_index,ket_orbital_index)*shell_coefficients(1)
                    elemnt6 = elemnt6 + dint6(orbital_index,ket_orbital_index)*shell_coefficients(1)
                    elemnt7 = elemnt7 + dint7(orbital_index,ket_orbital_index)*shell_coefficients(1)
                  end if
                end if
              end if
            end if
            do 9 j = 1,nvec
              loc = (j-1)*csf_count
              contri1 = csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)*elemnt1
              contri2 = csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)*elemnt2
              contri3 = csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)*elemnt3
              contri4 = csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)*elemnt4
              contri5 = csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)*elemnt5
              contri6 = csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)*elemnt6
              contri7 = csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)*elemnt7
              if (ir/=csf_index) then
                contri1 = 2.0_real64 * contri1
                contri2 = 2.0_real64 * contri2
                contri3 = 2.0_real64 * contri3
                contri4 = 2.0_real64 * contri4
                contri5 = 2.0_real64 * contri5
                contri6 = 2.0_real64 * contri6
                contri7 = 2.0_real64 * contri7
              end if
              dens1(j) = dens1(j) + contri1
              dens2(j) = dens2(j) + contri2
              dens3(j) = dens3(j) + contri3
              dens4(j) = dens4(j) + contri4
              dens5(j) = dens5(j) + contri5
              dens6(j) = dens6(j) + contri6
              dens7(j) = dens7(j) + contri7
    9       continue
          end if
   12   continue
   13 continue

       if (doit==1) write(50) -1
!
! Empty the buffer and close file
      if (doit==1) close(50)
!
!   Deallocate storage for the arrays in BUFFER
      call alcbuf (3)
      return
      end subroutine densnew
