      module read1_mem_i
      interface
!   Written by  G. Gaigalas               Vilnius, September 2021
      subroutine read1_mem(nfile, ncontr_tot2, lab, ncontr)
      use vast_kind_param, only: long
      integer(long), intent(in) :: ncontr_tot2
      integer, intent(in)       :: nfile
      integer, intent(out)      :: lab
      integer, intent(out)      :: ncontr
      end subroutine
      end interface
      end module
