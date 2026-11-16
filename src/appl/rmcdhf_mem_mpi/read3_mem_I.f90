      module read3_mem_i
      interface
!   Written by  G. Gaigalas               Vilnius, September 2021
      subroutine read3_mem(nfile,ncontr_tot,ncontr)
      use vast_kind_param, only: long
      integer(long), intent(in) :: ncontr_tot
      integer, intent(in)       :: nfile
      integer, intent(out)      :: ncontr
      end subroutine
      end interface
      end module
