      module setmcp_mem_i
      implicit none (type, external)
      interface
!   Written by  G. Gaigalas               Vilnius, September 2021
      subroutine setmcp_mem(nblock,nfile)
      implicit none (type, external)
      integer, intent(in) :: nblock
      integer, intent(in) :: nfile
      end subroutine setmcp_mem
      end interface
      end module setmcp_mem_i
