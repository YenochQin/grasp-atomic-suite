      module mcp_coefficient_storage
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer :: ncoeff, nintg
      integer, pointer :: jann(:), jbnn(:), intgrl(:), intptr(:)
      real(real64), pointer :: cnn(:)
      end module mcp_coefficient_storage
