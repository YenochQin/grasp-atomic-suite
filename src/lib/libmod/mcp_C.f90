      module mcp_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer :: max_coulomb_multipole_rank
      logical :: diagonal_only, use_first_order_restriction
      end module mcp_c
