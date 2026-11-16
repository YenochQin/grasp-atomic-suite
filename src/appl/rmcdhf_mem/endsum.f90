!***********************************************************************
!                                                                      *
      subroutine endsum
!                                                                      *
!   Generates the last part of  rscf92.sum  (on stream 24).            *
!                                                                      *
!   Call(s) to: [LIB92]: ENGOUT, RINT.                                 *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use eigv_c
      use orb_c
      use scf_c
!gg      USE syma_C
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rint_i
      use engoutgg_i
      use csfwgt_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, mode
      real(double) :: wa, wb, wc, wd, we
!-----------------------------------------------
!
!   Write out the orbital properties
!
      write (24, 301)
      do i = 1, nw
         write (24, 302) np(i),nh(i),e(i),pz(i),gama(i),pf(2,i),qf(2,i), &
            scnsty(i), mf(i)
      end do
!
      write (24, 303)
      do i = 1, nw
         wa = rint(i,i,-1)
         wb = rint(i,i,1)
         wc = rint(i,i,2)
         wd = rint (i,i, 4)
         we = 0.d0
         if (nh(i) /= 's ' .and. nh(i) /= 'p-') then
            we = rint(i,i,-3)
         end if
         write (24,304) np(i),nh(i),we,wa,wb,wc,wd, ucf(i)
      end do
!
      if (ncmin /= 0) then
         mode = 0
         call engoutgg (eval,iccmin,ncmin,mode)
!GG         CALL ENGOUT (EVAL, IATJPO, IASPAR, ICCMIN, NCMIN, MODE)
         call csfwgt (.false.)
      endif
!
      close(24)
!
      return
!
  301 format(/,'Radial wavefunction summary:'/,/,67x,'Self'/,'Subshell',6x,'e',&
         13x,'p0',5x,'gamma',5x,'P(2)',7x,'Q(2)',3x,'Consistency',' MTP'/)
  302 format(1x,i2,a2,1x,1p,d17.10,1p,d11.3,0p,f6.2,1p,3(d11.3),i5)
  303 format (/18x,'-3',14x,'-1',29x,'2',14x,'4',5x,'Generalised'     &
              /'Subshell',4x,'<  r  >',8x,'<  r  >',8x,'<  r  >',8x,  &
               '<  r  >',8x,'<  r  >',6x,'occupation'/)
  304 format (1x,1i2,1a2,1x,1p,6d15.5)
      return
!
      end subroutine endsum
