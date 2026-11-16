      module itjj2_i
      implicit none (type, external)
      interface
!
      integer function itjj2(ik,id,kg1,bk,bd,ibt,bt,itp,itg,iq)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)                :: kg1, iq
      integer,      intent(out)               :: itp, itg
      integer,      intent(in) :: ik(7), id(7)
      integer,      intent(out) :: ibt(7)
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out) :: bt(3)
      end function
      end interface
      end module
