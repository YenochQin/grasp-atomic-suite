!***********************************************************************
!                                                                      *
      integer function indtpi (itype, i)
!                                                                      *
!   Written by Farid A Parpia             Last revision: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
      use bilst_c, only: double, ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
          valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
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
