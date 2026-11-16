      module bess_c
      use vast_kind_param, only:  double
      use parameter_def,   only:  nnnp
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double), dimension(2) :: wij
      real(double), dimension(2,2,nnnp) :: bessj, bessn
      real(double), dimension(nnnp,3) :: bj
      real(double), dimension(nnnp) :: tc, td
      end module bess_c
