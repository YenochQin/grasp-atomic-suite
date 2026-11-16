!***********************************************************************
!                                                                      *
      subroutine densnew_mpi(doit,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
!                                                                      *
!   MPI version of DENSNEW. Each rank handles a strided subset of the  *
!   outer IC loop and can optionally write a per-rank angular-data      *
!   fragment that the root process later merges into name.IOB.          *
!                                                                      *
!***********************************************************************
      use vast_kind_param,  only: double
      use parameter_def,    only: orbital_index_radix, max_orbitals, max_radial_points
      use debug_c
      use decide_c
      use def_c
      use eigv_c
      use foparm_c
      use grid_c
      use jlabl_c
      use npar_c
      use orb_c
      use blk_c
      use prnt_c
      use teilst_c
      use buffer_c
      use ris_c
      use syma_c
      use mpi_c,            only: mpi_rank, mpi_process_count
      use prnt_c,           only : nvec
      use alcbuf_i
      use convrt_i
      use ispar_i
      use itjpo_i
      use onescalar_i
      implicit none
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint1, dint2, &
                                                        dint3, dint4, &
                                                        dint5, dint6, &
                                                        dint7
      integer, intent(in) :: doit
      integer, parameter :: key = orbital_index_radix
      real(double), dimension(max_orbitals) :: tshell
      real(double), dimension(max_orbitals) :: tshell_s
      real(double) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5,     &
                      elemnt6, elemnt7, contri1, contri2, contri3,     &
                      contri4, contri5, contri6, contri7, evpair,      &
                      pair_scale
      logical :: has_contr
      character :: cnum*11, ck*2
      integer, dimension(max_orbitals) :: ia_s
      integer :: ka, iopar, incor, ic, lcnum, itjpoc, ir, ia, ib, i, j
      integer :: loc, ncontr, lab, occ_ic, occ_ir, ndiff

      ka = 0
      iopar = 1
      incor = 1

      call alcbuf (1)

      do ic = mpi_rank + 1, csf_count, mpi_process_count
         if (mpi_rank == 0 .and. mod(ic,100) .eq. 1) then
            call convrt (ic,cnum,lcnum)
            print *, 'Column '//cnum(1:lcnum)//' complete on root;'
         endif

         itjpoc = itjpo (ic)

         do ir = ic,csf_count
            if (itjpo(ir) .eq. itjpoc) then
               if (ispar(ic) .ne. ispar(ir)) cycle
               ndiff = 0
               do i = 1,orbital_count
                  occ_ic = csf_subshell_occupations(i,ic)
                  occ_ir = csf_subshell_occupations(i,ir)
                  if (occ_ic .eq. occ_ir) cycle
                  if (iabs(occ_ic - occ_ir) .gt. 1) then
                     ndiff = 3
                     exit
                  endif
                  ndiff = ndiff + 1
                  if (ndiff .gt. 2) exit
               end do
               if (ndiff .ne. 0 .and. ndiff .ne. 2) cycle
               if (ndiff .eq. 0 .and. ic .ne. ir) cycle

               elemnt1 = 0.0d00
               elemnt2 = 0.0d00
               elemnt3 = 0.0d00
               elemnt4 = 0.0d00
               elemnt5 = 0.0d00
               elemnt6 = 0.0d00
               elemnt7 = 0.0d00
               has_contr = .false.

               call onescalar(ic,ir,ia,ib,tshell)
               if (ia .ne. 0) then
                  if (ia .eq. ib) then
                     ncontr = 0
                     do ia = 1,orbital_count
                        if (abs (tshell(ia)) .gt. cutoff) then
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
                        endif
                     end do
                     if (doit.eq.1) write(50) ic,ir,ncontr
                     if (doit.eq.1) then
                        do i = 1,ncontr
                           lab = ia_s(i)*(key + 1)
                           write(50) tshell_s(i),lab
                        end do
                     endif
                  else
                     if (abs (tshell(1)) .gt. cutoff) then
                        if (orbital_kappa(ia).eq.orbital_kappa(ib)) then
                           if (doit.eq.1) write(50) ic,ir,1
                           if (doit.eq.1) then
                              lab = ia*key + ib
                              write(50) tshell(1),lab
                           endif
                           elemnt1 = elemnt1 + dint1(ia,ib)*tshell(1)
                           elemnt2 = elemnt2 + dint2(ia,ib)*tshell(1)
                           elemnt3 = elemnt3 + dint3(ia,ib)*tshell(1)
                           elemnt4 = elemnt4 + dint4(ia,ib)*tshell(1)
                           elemnt5 = elemnt5 + dint5(ia,ib)*tshell(1)
                           elemnt6 = elemnt6 + dint6(ia,ib)*tshell(1)
                           elemnt7 = elemnt7 + dint7(ia,ib)*tshell(1)
                           has_contr = .true.
                        endif
                     endif
                  endif
               endif
               if (.not. has_contr) cycle
               if (abs(elemnt1) .le. cutoff .and. &
                   abs(elemnt2) .le. cutoff .and. &
                   abs(elemnt3) .le. cutoff .and. &
                   abs(elemnt4) .le. cutoff .and. &
                   abs(elemnt5) .le. cutoff .and. &
                   abs(elemnt6) .le. cutoff .and. &
                   abs(elemnt7) .le. cutoff) cycle
               if (ir .eq. ic) then
                  pair_scale = 1.0d00
               else
                  pair_scale = 2.0d00
               endif
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
            endif
         end do
      end do

      call alcbuf (3)
      return
      end subroutine densnew_mpi
