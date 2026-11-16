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
      use vast_kind_param,  only: double
      use parameter_def,    only: keyorb, nnnw, nnnp
      use blk_c
      use debug_c
      use decide_c
      use def_c
      use eigv_c
      use foparm_c
      use grid_c
      use jlabl_c
      use npar_c
      use orb_c
      use prnt_c
      use teilst_c
      use buffer_c
      use ris_c
      use syma_c
      use mpi_c,            only: myid, nprocs
      use prnt_c,           only : nvec
      use alcbuf_i
      use convrt_i
      use ispar_i
      use itjpo_i
      use onescalar_i
      implicit none
      integer, intent(in) :: doit, nrnuc
      real(double), dimension(nnnw,nnnw), intent(in) :: dint1, dint2, &
                                                        dint3, dint4, &
                                                        dint5, dint6, &
                                                        dint7
      real(double), dimension(nvec,nrnuc), intent(out) :: dens1vec
      real(double), dimension(nnnw,nnnw,nrnuc), intent(in) :: dint1vec
      integer, parameter :: key = keyorb
      real(double), dimension(nnnw) :: tshell, tshell_s
      real(double), dimension(nrnuc) :: contri1vec, elemnt1vec
      real(double) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5,     &
                      elemnt6, elemnt7, contri1, contri2, contri3,     &
                      contri4, contri5, contri6, contri7, evpair,      &
                      pair_scale
      logical :: has_contr
      character :: cnum*11, ck*2
      integer, dimension(nnnw) :: ia_s
      integer :: ka, iopar, incor, ic, lcnum, itjpoc, ir, ia, ib, i, j
      integer :: l, loc, ncontr, lab, occ_ic, occ_ir, ndiff

      dens1vec(:,:) = 0.0d00
      ka = 0
      iopar = 1
      incor = 1

      call alcbuf (1)

      do ic = myid + 1, ncf, nprocs
         if (myid == 0 .and. mod(ic,100) .eq. 1) then
            call convrt (ic,cnum,lcnum)
            print *, 'Column '//cnum(1:lcnum)//' complete on root;'
         endif

         itjpoc = itjpo (ic)
         do ir = ic,ncf
            if (itjpo(ir) .eq. itjpoc) then
               if (ispar(ic) .ne. ispar(ir)) cycle
               ndiff = 0
               do i = 1,nw
                  occ_ic = iqa(i,ic)
                  occ_ir = iqa(i,ir)
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
               elemnt1vec(:) = 0.0d00
               has_contr = .false.

               call onescalar(ic,ir,ia,ib,tshell)
               if (ia .ne. 0) then
                  if (ia .eq. ib) then
                     ncontr = 0
                     do ia = 1,nw
                        if (abs (tshell(ia)) .gt. cutoff) then
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
                        endif
                     end do
                  else
                     if (abs (tshell(1)) .gt. cutoff) then
                        if (nak(ia).eq.nak(ib)) then
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
                   abs(elemnt7) .le. cutoff .and. &
                   all(abs(elemnt1vec(:)) .le. cutoff)) cycle
               if (ir .eq. ic) then
                  pair_scale = 1.0d00
               else
                  pair_scale = 2.0d00
               endif
               do j = 1,nvec
                  loc = (j-1)*ncf
                  evpair = pair_scale * evec(ic+loc)*evec(ir+loc)
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
            endif
         end do
      end do

      call alcbuf (3)
      return
      end subroutine densnew_seltz_mpi
