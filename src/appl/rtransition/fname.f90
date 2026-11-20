!***********************************************************************
!                                                                      *
      subroutine fname(name)
      use calculation_state_files, only: read_state_pair
!                                                                      *
!   Determines the name of the initial and final states                *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=24), intent(out) :: name(2)
      call read_state_pair(name)

      return
      end subroutine fname
