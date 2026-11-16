      module izas1_i
      implicit none (type, external)
      interface
!
      integer function izas1(ib,qb,ik,qk)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in) :: ib, ik
      real(real64), intent(in) :: qb, qk
      end function
      end interface
      end module
