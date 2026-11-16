!***********************************************************************
!                                                                      *
      integer function ichkq1 (ja, jb)
!                                                                      *
!   This routine is to check the occupation condition for one electron *
!   operator.                                                          *
!                                                                      *
!   Call(s) to: [LIB92]: IQ.                                           *
!                                                                      *
!   Yu Zou                                Last revision: 8/16/00       *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:20   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use debug_c
      use orb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: ja
      integer :: jb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, i_iqa, i_iqb

!-----------------------------------------------
!
!
      ichkq1 = 0
      k = 0
      do i = 1, nw
         i_iqa = iq(i,ja)
         i_iqb = iq(i,jb)
         if (i_iqa == i_iqb) cycle
         k = k + 1
         if (k > 2) return
         if (iabs(i_iqa - i_iqb) <= 1) cycle
         return
      end do
      if (k==2 .or. k==0) ichkq1 = 1
      return
      end function ichkq1
