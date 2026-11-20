      module dsubrs_i
      implicit none (type, external)
      interface
      real(real64) function dsubrs (extended_optimal_level, bra_csf_index, ket_csf_index, symmetry_block_index)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      logical, intent(in) :: extended_optimal_level
      integer, intent(in) :: bra_csf_index
      integer, intent(in) :: ket_csf_index
      integer :: symmetry_block_index
      end function dsubrs
      end interface
      end module dsubrs_i
