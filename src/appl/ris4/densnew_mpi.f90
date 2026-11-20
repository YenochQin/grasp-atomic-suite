!***********************************************************************
!                                                                      *
      subroutine densnew_mpi(doit,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
!                                                                      *
!   MPI version of DENSNEW. Each rank handles a strided subset of the  *
!   outer IC loop and can optionally write a per-rank angular-data      *
!   fragment that the root process later merges into name.IOB.          *
!                                                                      *
!***********************************************************************
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
      use angular_symmetry_labels, only: jlbl, jlbp, jlbr
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count, orbital_kappa
      use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
          nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use teilst_c, only: ndtea, ntei, pindte, pvalte, firstt
      use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
          angular_coefficient_labels, angular_coefficients
      use ris_c, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use mpi_state,            only: mpi_rank, mpi_process_count
      use selected_level_indices,           only : nvec
      use angular_coefficient_buffer, only: alcbuf
      use integer_format, only: convrt
      use csf_quantum_queries, only: ispar, itjpo
      use one_electron_scalar_coefficients, only: onescalar
      implicit none (type, external)
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), &
                                                        dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                                        dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals), &
                                                        dint7(max_orbitals, max_orbitals)
      integer, intent(in) :: doit
      integer, parameter :: key = orbital_index_radix
      real(real64) :: tshell(max_orbitals)
      real(real64) :: tshell_s(max_orbitals)
      real(real64) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5,     &
                      elemnt6, elemnt7, contri1, contri2, contri3,     &
                      contri4, contri5, contri6, contri7, evpair,      &
                      pair_scale
      logical :: has_contr
      character :: cnum*11, ck*2
      integer :: ia_s(max_orbitals)
      integer :: ka, iopar, incor, ic, lcnum, itjpoc, ir, ia, ib, i, j
      integer :: loc, ncontr, lab, occ_ic, occ_ir, ndiff

      ka = 0
      iopar = 1
      incor = 1

      call alcbuf (1)

      do ic = mpi_rank + 1, csf_count, mpi_process_count
         if (mpi_rank == 0 .and. mod(ic,100) == 1) then
            call convrt (ic,cnum,lcnum)
            print *, "Column "//cnum(1:lcnum)//" complete on root;"
         end if

         itjpoc = itjpo (ic)

         do ir = ic,csf_count
            if (itjpo(ir) == itjpoc) then
               if (ispar(ic) /= ispar(ir)) cycle
               ndiff = 0
               do i = 1,orbital_count
                  occ_ic = csf_subshell_occupations(i,ic)
                  occ_ir = csf_subshell_occupations(i,ir)
                  if (occ_ic == occ_ir) cycle
                  if (iabs(occ_ic - occ_ir) > 1) then
                     ndiff = 3
                     exit
                  end if
                  ndiff = ndiff + 1
                  if (ndiff > 2) exit
               end do
               if (ndiff /= 0 .and. ndiff /= 2) cycle
               if (ndiff == 0 .and. ic /= ir) cycle

               elemnt1 = 0.0_real64
               elemnt2 = 0.0_real64
               elemnt3 = 0.0_real64
               elemnt4 = 0.0_real64
               elemnt5 = 0.0_real64
               elemnt6 = 0.0_real64
               elemnt7 = 0.0_real64
               has_contr = .false.

               call onescalar(ic,ir,ia,ib,tshell)
               if (ia /= 0) then
                  if (ia == ib) then
                     ncontr = 0
                     do ia = 1,orbital_count
                        if (abs (tshell(ia)) > cutoff) then
                           ncontr = ncontr + 1
                           tshell_s(ncontr) = tshell(ia)
                           ia_s(ncontr) = ia
                           elemnt1 = elemnt1 + dint1(ia,ia)*tshell(ia)
                           elemnt2 = elemnt2 + dint2(ia,ia)*tshell(ia)
                           elemnt3 = elemnt3 + dint3(ia,ia)*tshell(ia)
                           elemnt4 = elemnt4 + dint4(ia,ia)*tshell(ia)
                           elemnt5 = elemnt5 + dint5(ia,ia)*tshell(ia)
                           elemnt6 = elemnt6 + dint6(ia,ia)*tshell(ia)
                           elemnt7 = elemnt7 + dint7(ia,ia)*tshell(ia)
                           has_contr = .true.
                        end if
                     end do
                     if (doit==1) write(50) ic,ir,ncontr
                     if (doit==1) then
                        do i = 1,ncontr
                           lab = ia_s(i)*(key + 1)
                           write(50) tshell_s(i),lab
                        end do
                     end if
                  else
                     if (abs (tshell(1)) > cutoff) then
                        if (orbital_kappa(ia)==orbital_kappa(ib)) then
                           if (doit==1) write(50) ic,ir,1
                           if (doit==1) then
                              lab = ia*key + ib
                              write(50) tshell(1),lab
                           end if
                           elemnt1 = elemnt1 + dint1(ia,ib)*tshell(1)
                           elemnt2 = elemnt2 + dint2(ia,ib)*tshell(1)
                           elemnt3 = elemnt3 + dint3(ia,ib)*tshell(1)
                           elemnt4 = elemnt4 + dint4(ia,ib)*tshell(1)
                           elemnt5 = elemnt5 + dint5(ia,ib)*tshell(1)
                           elemnt6 = elemnt6 + dint6(ia,ib)*tshell(1)
                           elemnt7 = elemnt7 + dint7(ia,ib)*tshell(1)
                           has_contr = .true.
                        end if
                     end if
                  end if
               end if
               if (.not. has_contr) cycle
               if (abs(elemnt1) <= cutoff .and. &
                   abs(elemnt2) <= cutoff .and. &
                   abs(elemnt3) <= cutoff .and. &
                   abs(elemnt4) <= cutoff .and. &
                   abs(elemnt5) <= cutoff .and. &
                   abs(elemnt6) <= cutoff .and. &
                   abs(elemnt7) <= cutoff) cycle
               if (ir == ic) then
                  pair_scale = 1.0_real64
               else
                  pair_scale = 2.0_real64
               end if
               do j = 1,nvec
                  loc = (j-1)*csf_count
                  evpair = pair_scale * csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)
                  contri1 = evpair*elemnt1
                  contri2 = evpair*elemnt2
                  contri3 = evpair*elemnt3
                  contri4 = evpair*elemnt4
                  contri5 = evpair*elemnt5
                  contri6 = evpair*elemnt6
                  contri7 = evpair*elemnt7
                  dens1(j) = dens1(j) + contri1
                  dens2(j) = dens2(j) + contri2
                  dens3(j) = dens3(j) + contri3
                  dens4(j) = dens4(j) + contri4
                  dens5(j) = dens5(j) + contri5
                  dens6(j) = dens6(j) + contri6
                  dens7(j) = dens7(j) + contri7
               end do
            end if
         end do
      end do

      call alcbuf (3)
      return
      end subroutine densnew_mpi
