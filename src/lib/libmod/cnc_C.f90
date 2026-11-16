      module cnc_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      real(real64) :: cnc5c(2:5, 2:4)
      real(real64) :: cnc6c(6, 2:6)
      end module cnc_c
