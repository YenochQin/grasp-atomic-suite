!***********************************************************************
!                                                                      *
      subroutine densread_seltz_mpi(dint1,dint2,dint3,dint4,dint5,dint6,&
         dint7,dint1vec,dens1vec,nrnuc)
!                                                                      *
!   MPI version of DENSREAD_SELTZ. Each rank scans the full angular    *
!   file but only accumulates records whose IC belongs to its strided   *
!   subset.                                                             *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: orbital_index_radix, max_orbitals, max_radial_points
      use debug_output_settings,          only: cutoff
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use ris_c, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
      use orbital_csf_basis, only: csf_count, max_orbitals
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use mpi_state,            only: mpi_rank, mpi_process_count
      implicit none (type, external)
      integer, intent(in) :: nrnuc
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals),  &
                                          dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                             dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals),  &
                                          dint7(max_orbitals, max_orbitals)
      real(real64), intent(out) :: dens1vec(nvec, nrnuc)
      real(real64), intent(in) :: dint1vec(max_orbitals, max_orbitals, nrnuc)
      integer, parameter :: key = orbital_index_radix
      real(real64) :: tshell_r(max_orbitals)
      real(real64) :: contri1vec(nrnuc), elemnt1vec(nrnuc)
      real(real64) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5
      real(real64) :: elemnt6, elemnt7
      real(real64) :: contri1, contri2, contri3, contri4, contri5
      real(real64) :: contri6, contri7
      real(real64) :: evpair, pair_scale
      integer :: ios, ia, ib, ic, ir, i, j, l, loc, lab, ncount
      logical :: ownrec

      dens1vec(:,:) = 0.0_real64
      rewind (50)
   16 read (50,iostat = ios) ic,ir,ncount

      if (ios == 0) then
         ownrec = mod(ic-1, mpi_process_count) == mpi_rank
         elemnt1 = 0.0_real64
         elemnt2 = 0.0_real64
         elemnt3 = 0.0_real64
         elemnt4 = 0.0_real64
         elemnt5 = 0.0_real64
         elemnt6 = 0.0_real64
         elemnt7 = 0.0_real64
         elemnt1vec(:) = 0.0_real64

         do i= 1,ncount
            read(50) tshell_r(i),lab
            if (.not. ownrec) cycle
            ia = lab/key
            ib = mod(lab,key)
            elemnt1 = elemnt1 + dint1(ia,ib)*tshell_r(i)
            do l = 2,nrnuc
               elemnt1vec(l) = elemnt1vec(l) + dint1vec(ia,ib,l)*tshell_r(i)
            end do
            elemnt2 = elemnt2 + dint2(ia,ib)*tshell_r(i)
            elemnt3 = elemnt3 + dint3(ia,ib)*tshell_r(i)
            elemnt4 = elemnt4 + dint4(ia,ib)*tshell_r(i)
            elemnt5 = elemnt5 + dint5(ia,ib)*tshell_r(i)
            elemnt6 = elemnt6 + dint6(ia,ib)*tshell_r(i)
            elemnt7 = elemnt7 + dint7(ia,ib)*tshell_r(i)
         end do

         if (ownrec) then
            if (abs(elemnt1) <= cutoff .and. &
                abs(elemnt2) <= cutoff .and. &
                abs(elemnt3) <= cutoff .and. &
                abs(elemnt4) <= cutoff .and. &
                abs(elemnt5) <= cutoff .and. &
                abs(elemnt6) <= cutoff .and. &
                abs(elemnt7) <= cutoff .and. &
                all(abs(elemnt1vec(:)) <= cutoff)) goto 16
            if (ir == ic) then
               pair_scale = 1.0_real64
            else
               pair_scale = 2.0_real64
            end if
            do j = 1,nvec
               loc = (j-1)*csf_count
               evpair = pair_scale * csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)
               contri1 = evpair*elemnt1
               contri1vec(:) = evpair*elemnt1vec(:)
               contri2 = evpair*elemnt2
               contri3 = evpair*elemnt3
               contri4 = evpair*elemnt4
               contri5 = evpair*elemnt5
               contri6 = evpair*elemnt6
               contri7 = evpair*elemnt7
               dens1(j) = dens1(j) + contri1
               do l = 2,nrnuc
                  dens1vec(j,l) = dens1vec(j,l) + contri1vec(l)
               end do
               dens2(j) = dens2(j) + contri2
               dens3(j) = dens3(j) + contri3
               dens4(j) = dens4(j) + contri4
               dens5(j) = dens5(j) + contri5
               dens6(j) = dens6(j) + contri6
               dens7(j) = dens7(j) + contri7
            end do
         end if

         goto 16
      end if

      return
      end subroutine densread_seltz_mpi
