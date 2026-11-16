      module dumx_c
      use cons_c
!...Created by Pacific-Sierra Research 77to90  4.3E  10:50:24   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      use parameter_def,   only:  nnnw
      integer, dimension(nnnw) :: jlis, jc1s, jc2s
!$OMP THREADPRIVATE(JLIS, JC1S, JC2S)
      end module dumx_c
