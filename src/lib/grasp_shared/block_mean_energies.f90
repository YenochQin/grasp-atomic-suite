      module block_mean_energies
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      real(real64), pointer :: eavblk(:)
      end module block_mean_energies
