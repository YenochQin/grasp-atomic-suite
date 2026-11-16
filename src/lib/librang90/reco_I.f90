      module reco_i
      interface
!
      subroutine reco(ja1,ja2,ja3,ja4,ka,iat)
      use vast_kind_param, only:  double
      integer, intent(in)  :: ja1, ja2, ja3, ja4, ka
      integer, intent(out) :: iat
      end subroutine
      end interface
      end module
