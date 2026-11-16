      module csfm_i
      implicit none (type, external)
      interface
      subroutine csfm (asfa,asfb,lev1,lev2)
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)
        real(real64), intent(out) :: asfa, asfb
        integer, intent(in) :: lev1, lev2
      end subroutine csfm
      end interface
      end module csfm_i
