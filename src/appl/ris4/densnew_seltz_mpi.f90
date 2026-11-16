!***********************************************************************
!                                                                      *
      subroutine densnew_seltz_mpi(doit,dint1,dint2,dint3,dint4,dint5, &
         dint6,dint7,dint1vec,dens1vec,nrnuc)
!                                                                      *
!   MPI version of DENSNEW_SELTZ. Each rank handles a strided subset   *
!   of the outer IC loop. In MPI mode angular coefficients are         *
!   evaluated on the fly and not saved.                                *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: orbital_index_radix, max_orbitals, max_radial_points
      use blk_c
      use debug_c
      use decide_c
      use def_c
      use eigv_c
      use foparm_c, only: iccut
      use grid_c
      use jlabl_c, only: jlbl, jlbp, jlbr
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use orb_c
      use prnt_c
      use teilst_c
      use buffer_c
      use ris_c
      use syma_c
      use mpi_c,            only: mpi_rank, mpi_process_count
      use prnt_c,           only : nvec
      use alcbuf_i, only: alcbuf
      use convrt_i, only: convrt
      use ispar_i, only: ispar
      use itjpo_i, only: itjpo
      use onescalar_i, only: onescalar
      implicit none (type, external)
      integer, intent(in) :: doit, nrnuc
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), &
                                                        dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                                        dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals), &
                                                        dint7(max_orbitals, max_orbitals)
      real(real64), intent(out) :: dens1vec(nvec, nrnuc)
      real(real64), intent(in) :: dint1vec(max_orbitals, max_orbitals, nrnuc)
      integer, parameter :: key = orbital_index_radix
      real(real64) :: tshell(max_orbitals), tshell_s(max_orbitals)
      real(real64) :: contri1vec(nrnuc), elemnt1vec(nrnuc)
      real(real64) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5,     &
                      elemnt6, elemnt7, contri1, contri2, contri3,     &
                      contri4, contri5, contri6, contri7, evpair,      &
                      pair_scale
      logical :: has_contr
      character :: cnum*11, ck*2
      integer :: ia_s(max_orbitals)
      integer :: ka, iopar, incor, ic, lcnum, itjpoc, ir, ia, ib, i, j
      integer :: l, loc, ncontr, lab, occ_ic, occ_ir, ndiff

      dens1vec(:,:) = 0.0_real64
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
               elemnt1vec(:) = 0.0_real64
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
                           do l = 2,nrnuc
                              elemnt1vec(l)=elemnt1vec(l)+dint1vec(ia,ia,l)*tshell(ia)
                           end do
                           elemnt2 = elemnt2 + dint2(ia,ia)*tshell(ia)
                           elemnt3 = elemnt3 + dint3(ia,ia)*tshell(ia)
                           elemnt4 = elemnt4 + dint4(ia,ia)*tshell(ia)
                           elemnt5 = elemnt5 + dint5(ia,ia)*tshell(ia)
                           elemnt6 = elemnt6 + dint6(ia,ia)*tshell(ia)
                           elemnt7 = elemnt7 + dint7(ia,ia)*tshell(ia)
                           has_contr = .true.
                        end if
                     end do
                  else
                     if (abs (tshell(1)) > cutoff) then
                        if (orbital_kappa(ia)==orbital_kappa(ib)) then
                           elemnt1 = elemnt1 + dint1(ia,ib)*tshell(1)
                           do l = 2,nrnuc
                              elemnt1vec(l)=elemnt1vec(l)+dint1vec(ia,ib,l)*tshell(1)
                           end do
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
                   abs(elemnt7) <= cutoff .and. &
                   all(abs(elemnt1vec(:)) <= cutoff)) cycle
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
         end do
      end do

      call alcbuf (3)
      return
      end subroutine densnew_seltz_mpi
