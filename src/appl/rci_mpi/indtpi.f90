!***********************************************************************
!                                                                      *
      integer function indtpi (itype, i)
!                                                                      *
!   Written by Farid A Parpia             Last revision: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      use bilst_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: itype
      integer  :: i
!-----------------------------------------------
!
      select case (itype)
      case (1)
         indtpi = indtp1(i)
      case (2)
         indtpi = indtp2(i)
      case (3)
         indtpi = indtp3(i)
      case (4)
         indtpi = indtp4(i)
      case (5)
         indtpi = indtp5(i)
      case (6)
         indtpi = indtp6(i)
      end select
!
      return
      end function indtpi
