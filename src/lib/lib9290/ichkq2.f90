!***********************************************************************
!                                                                      *
      integer function ichkq2 (ja, jb)
!                                                                      *
!   This routine is to check the occupation condition for two electron *
!   operator.                                                          *
!                                                                      *
!   Call(s) to: [LIB92]: IQ.                                           *
!                                                                      *
!   Yu Zou                                Last revision: 8/21/00       *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_c
      use orb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ja
      integer, intent(in) :: jb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, i_iqa, i_iqb
!-----------------------------------------------
!
!
      ichkq2 = 0
      k = 0
      do i = 1, orbital_count
         i_iqa = iq(i,ja)
         i_iqb = iq(i,jb)
         if (i_iqa == i_iqb) cycle
         k = k + 1
         if (k > 4) return
         if (iabs(i_iqa - i_iqb) <= 2) cycle
         return
      end do
      ichkq2 = 1
      return
      end function ichkq2
