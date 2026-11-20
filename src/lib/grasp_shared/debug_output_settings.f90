! Stores switches for writes from the angular, radial, and property routines.
      module debug_output_settings
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer :: ibug1, ibug2, ibug3, ibug4, ibug5, ibug6
      logical :: ldbpa(5)
      logical :: ldbpg(5)
      logical :: ldbpr(30)
      real(real64) :: cutoff  ! Threshold for the matrix-element trace in BIOSCL.
      end module debug_output_settings
