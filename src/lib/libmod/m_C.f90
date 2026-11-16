!
!***********************************************************************
!                                                                      *
      module m_c
      use vast_kind_param, only:  double
      use parameter_def,   only:  max_orbitals
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(max_orbitals) :: nq1, nq2
      integer, dimension(max_orbitals) :: jjc1, jjc2
      integer, dimension(3,max_orbitals) :: jjq1, jjq2
      integer, dimension(max_orbitals) :: jlist, klist
      integer :: npeel, ncore
!$OMP THREADPRIVATE(NQ1, NQ2, JJC1, JJC2, JJQ1, JJQ2, JLIST, KLIST, NPEEL, NCORE)
      end module m_c
