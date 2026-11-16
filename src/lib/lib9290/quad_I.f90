      module quad_i
      implicit none (type, external)
      interface
      subroutine quad (result)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(out) :: result
!VAST.../GRID/ R(IN), RP(IN)
!VAST.../NCC/ C1(IN), C2(IN), C3(IN), C4(IN)
!VAST.../TATB/ TA(INOUT), MTP(IN)
      end subroutine
      end interface
      end module
