!***********************************************************************
!                                                                      *
!JE   SUBROUTINE DENSNEW(DOIT,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
       subroutine densnew_seltz(doit,dint1,dint2,dint3,                &
                          dint4,dint5,dint6,dint7,                     &
                          dint1vec,dens1vec,nrnuc)
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
      use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
          nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use atomic_calculation_state, only: max_radial_points
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use first_order_csf_cutoffs, only: iccut
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use angular_symmetry_labels, only: jlbl, jlbp, jlbr
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count, orbital_kappa
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
      integer :: l, loc, ncontr, lab
      integer, intent(in) :: nrnuc
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), &
                                                        dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                                        dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals), &
                                                        dint7(max_orbitals, max_orbitals)
      real(real64), intent(out)     :: dens1vec(nvec, nrnuc) !  JE ADD
      real(real64), intent(in) :: dint1vec(max_orbitals, max_orbitals, nrnuc) !  JE ADD
      integer, intent(in) :: doit
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: tshell(max_orbitals), tshell_s(max_orbitals)
      real(real64) :: contri1vec(nrnuc), elemnt1vec(nrnuc)
      real(real64), pointer :: emt1(:), emt2(:), emt3(:), emt4(:),   &
                                             emt5(:), emt6(:)
      real(real64) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5,     &
                      elemnt6, elemnt7, contri1, contri2, contri3,     &
                      contri4, contri5, contri6, contri7
      logical :: vsh, nucde, smssh, yes
      character :: cnum*11, ck*2
      integer :: ia_s(max_orbitals)
      integer :: ka, iopar, incor, ic, lcnum, itjpoc, ir, ia, ib, i, ii, j
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
      dens1vec(:,:) = 0.0_real64                                     ! JE ADD
!
!   Set the rank (zero) and parity (even) for the one-particle
!   coefficients
!
      ka = 0
      iopar = 1
      incor = 1
!
!   Allocate storage for the arrays in BUFFER
!
      call alcbuf (1)
!
!   Sweep through the Hamiltonian matrix to determine the
!   sms parameter
!
         do 13 ic = 1,csf_count
!
!   Output IC on the screen to show how far the calculation has preceede
!
        call convrt (ic,cnum,lcnum)
        if (mod(ic,100)==0) then
          print *, "Column "//cnum(1:lcnum)//" complete;"
        end if
!
        itjpoc = itjpo (ic)

        do 12 ir = ic,csf_count
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
            elemnt1vec(:) = 0.0_real64                                      ! JE ADD
!
!   Call the MCT package to compute T coefficients
!
           call onescalar(ic,ir,ia,ib,tshell)
!GG            CALL TNSRJJ (KA,IOPAR,IC,IR,IA,IB,TSHELL)
            if (ia /= 0) then
              if (ia == ib) then
                    ncontr = 0
                do 8 ia = 1,orbital_count
                  if (abs (tshell(ia)) > cutoff) then
                    ncontr = ncontr + 1
                    tshell_s(ncontr) = tshell(ia)
                    ia_s(ncontr) = ia
                    elemnt1 = elemnt1 + dint1(ia,ia)*tshell(ia)
                    do 21 l = 2,nrnuc                                    ! JE ADD
                      elemnt1vec(l) = elemnt1vec(l)+dint1vec(ia,ia,l)*& ! JE ADD
                                      tshell(ia)                        ! JE ADD
   21               continue                                            ! JE ADD
                    elemnt2 = elemnt2 + dint2(ia,ia)*tshell(ia)
                    elemnt3 = elemnt3 + dint3(ia,ia)*tshell(ia)
                    elemnt4 = elemnt4 + dint4(ia,ia)*tshell(ia)
                    elemnt5 = elemnt5 + dint5(ia,ia)*tshell(ia)
                    elemnt6 = elemnt6 + dint6(ia,ia)*tshell(ia)
                    elemnt7 = elemnt7 + dint7(ia,ia)*tshell(ia)
                  end if
    8           continue
                if (doit==1) write(50) ic,ir,ncontr
                do i = 1,ncontr
                  lab = ia_s(i)*(key + 1)
                  if (doit==1) write(50) tshell_s(i),lab
                end do
              else
                if (abs (tshell(1)) > cutoff) then
                  if (orbital_kappa(ia)==orbital_kappa(ib)) then
                    if (doit==1) write(50) ic,ir,1
                    lab = ia*key + ib
                    if (doit==1) write(50) tshell(1),lab
                    elemnt1 = elemnt1 + dint1(ia,ib)*tshell(1)
                    do 23 l = 2,nrnuc                                   ! JE ADD
                      elemnt1vec(l) = elemnt1vec(l)+dint1vec(ia,ib,l)*& ! JE ADD
                                      tshell(1)                         ! JE ADD
   23               continue                                            ! JE ADD
                    elemnt2 = elemnt2 + dint2(ia,ib)*tshell(1)
                    elemnt3 = elemnt3 + dint3(ia,ib)*tshell(1)
                    elemnt4 = elemnt4 + dint4(ia,ib)*tshell(1)
                    elemnt5 = elemnt5 + dint5(ia,ib)*tshell(1)
                    elemnt6 = elemnt6 + dint6(ia,ib)*tshell(1)
                    elemnt7 = elemnt7 + dint7(ia,ib)*tshell(1)
                  end if
                end if
              end if
            end if
            do 9 j = 1,nvec
              loc = (j-1)*csf_count
              contri1 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt1
              contri1vec(:) = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt1vec(:)   ! JE ADD
              contri2 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt2
              contri3 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt3
              contri4 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt4
              contri5 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt5
              contri6 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt6
              contri7 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt7
              if (ir/=ic) then
                contri1 = 2.0_real64 * contri1
                contri1vec(:) = 2.0_real64 * contri1vec(:)                  ! JE ADD
                contri2 = 2.0_real64 * contri2
                contri3 = 2.0_real64 * contri3
                contri4 = 2.0_real64 * contri4
                contri5 = 2.0_real64 * contri5
                contri6 = 2.0_real64 * contri6
                contri7 = 2.0_real64 * contri7
              end if
              dens1(j) = dens1(j) + contri1
              do 22 l = 2,nrnuc                                         ! JE ADD
                dens1vec(j,l) = dens1vec(j,l) + contri1vec(l)           ! JE ADD
   22         continue                                                  ! JE ADD
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
      end subroutine densnew_seltz
