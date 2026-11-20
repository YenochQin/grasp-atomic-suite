      module ielsum_i
      implicit none (type, external)
      interface
      integer function ielsum (ivec, nelmnt)
      implicit none (type, external)
      integer, intent(in) :: ivec(*)
      integer, intent(in) :: nelmnt
      end function ielsum
      end interface
      end module ielsum_i
