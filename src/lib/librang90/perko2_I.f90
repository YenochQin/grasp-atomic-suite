      module perko2_i
      interface
      subroutine perko2(ja1,ja2,ja3,ja4,i)
      use vast_kind_param, only:  double
      use trk_c
      integer, intent(in) :: ja1, ja2, ja3, ja4, i
      end subroutine
      end interface
      end module
