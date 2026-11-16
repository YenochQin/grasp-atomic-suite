!
!***********************************************************************
!                                                                      *
      module debug_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: ibug1, ibug2, ibug3, ibug4, ibug5, ibug6
      logical :: ldbpa(5)
      logical :: ldbpg(5)
      logical :: ldbpr(30)
      real(real64) :: cutoff  ! used by bioscl
      end module debug_c
