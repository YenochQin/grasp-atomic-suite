!
!***********************************************************************
!                                                                      *
      module m_c
      use vast_kind_param, only:  double
      use parameter_def,   only:  nnnw
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(nnnw) :: nq1, nq2
      integer, dimension(nnnw) :: jjc1, jjc2
      integer, dimension(3,nnnw) :: jjq1, jjq2
      integer, dimension(nnnw) :: jlist, klist
      integer :: npeel, ncore
!$OMP THREADPRIVATE(NQ1, NQ2, JJC1, JJC2, JJQ1, JJQ2, JLIST, KLIST, NPEEL, NCORE)
      end module m_c
