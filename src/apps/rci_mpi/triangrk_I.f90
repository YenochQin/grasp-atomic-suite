      module triangrk_i
      implicit none (type, external)
      interface
      logical function triangrk (la, k, lb)
      implicit none (type, external)
      integer, intent(in) :: la
      integer, intent(in) :: k
      integer, intent(in) :: lb
      end function triangrk
      end interface
      end module triangrk_i
