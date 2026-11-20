!***********************************************************************

      subroutine setsum(name)
      use calculation_file_units, only: summary_output_unit
!
!   Open the  .csum  file on stream 24.
!   Xinghong He                                          10 Jun 1998
!
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_open, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, open_status
!-----------------------------------------------------------------------

      k = index(name," ")

      call openfl (summary_output_unit, name(1:k-1)//".csum", "FORMATTED", "UNKNOWN", open_status)

      return
      end subroutine setsum
