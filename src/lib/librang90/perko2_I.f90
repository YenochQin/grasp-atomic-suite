      module perko2_i
      implicit none (type, external)
      interface
      subroutine perko2(ja1,ja2,ja3,ja4,i)
      use trk_c
      implicit none (type, external)
      integer, intent(in) :: ja1, ja2, ja3, ja4, i
      end subroutine
      end interface
      end module
