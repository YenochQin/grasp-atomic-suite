      module read2_mem_i
      interface
!   Written by  G. Gaigalas               Vilnius, September 2021
      subroutine read2_mem(nfile,ncontr_tot,ncontr,iclmndum)
      use vast_kind_param, only: long
      integer(long), intent(in) :: ncontr_tot
      integer, intent(in)       :: nfile
      integer, intent(out)      :: ncontr
      integer, intent(out)      :: iclmndum
      end subroutine
      end interface
      end module
