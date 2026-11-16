      module read1_mem_i
      implicit none (type, external)
      interface
!   Written by  G. Gaigalas               Vilnius, September 2021
      subroutine read1_mem(nfile, ncontr_tot2, lab, ncontr)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none (type, external)
      integer(int64), intent(in) :: ncontr_tot2
      integer, intent(in)       :: nfile
      integer, intent(out)      :: lab
      integer, intent(out)      :: ncontr
      end subroutine read1_mem
      end interface
      end module read1_mem_i
