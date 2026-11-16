!***********************************************************************
!
      subroutine getald
!
!   Interactively determines the data governing AL problem.
!
!   Quantities to be determined:
!     ncmin, wt(1:ncf), ucf(1:nw), nscf, nsic, orthst
!
!   Call(s) to: [LIB92]: ALLOC, IQ.
!
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992
!   Block version by Xinghong He          Last revision: 13 Jul 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use def_c
      use fixd_c
      use orb_c
      use orthct_c
      use scf_c
      use iounit_c
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getaldwt_i
      use iq_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqadum, j, i
      real(double) :: sum
      logical :: getyn, yes
!-----------------------------------------------
!
      write (istde, *) '(E)AL type calculation; H(DC) will not be ', &
         'diagonalised;'
      write (istde, *) 'getald ...'
      write (istde, *) 'ncf=', ncf

      call alloc (wt, ncf, 'WT', 'GETALD')

      call getaldwt (ncf, wt)

      do j = 1, nw
         sum = 0.d0
         do i = 1, ncf
            sum = sum + wt(i)*dble(iq(j,i))
         end do
         ucf(j) = sum
      end do

      ncmin = 0
      nscf = 12
      nsic = 2 + (nw - nfix)/4
      orthst = .false.

      return
      end subroutine getald
