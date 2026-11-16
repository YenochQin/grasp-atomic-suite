      module eigvec1_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      real(real64), pointer :: evec1(:)
      end module eigvec1_c
