      module recop00_i
      interface
!
      subroutine recop00(ns,ja1,ja2,ka,iat)
      use vast_kind_param, only:  double
      integer, intent(in)  :: ns, ja1, ja2, ka
      integer, intent(out) :: iat
      end subroutine
      end interface
      end module
