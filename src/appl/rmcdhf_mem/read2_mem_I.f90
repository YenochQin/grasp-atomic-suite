      module read2_mem_i
      implicit none (type, external)
      interface
!   Written by  G. Gaigalas               Vilnius, September 2021
      subroutine read2_mem(nfile,ncontr_tot,ncontr,iclmndum)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none (type, external)
      integer(int64), intent(in) :: ncontr_tot
      integer, intent(in)       :: nfile
      integer, intent(out)      :: ncontr
      integer, intent(out)      :: iclmndum
      end subroutine read2_mem
      end interface
      end module read2_mem_i
