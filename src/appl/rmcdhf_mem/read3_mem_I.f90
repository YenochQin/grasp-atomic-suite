      module read3_mem_i
      implicit none (type, external)
      interface
!   Written by  G. Gaigalas               Vilnius, September 2021
      subroutine read3_mem(nfile,ncontr_tot,ncontr)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none (type, external)
      integer(int64), intent(in) :: ncontr_tot
      integer, intent(in)       :: nfile
      integer, intent(out)      :: ncontr
      end subroutine
      end interface
      end module
