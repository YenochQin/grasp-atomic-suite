      module breid_i
      implicit none (type, external)
      interface
      subroutine breid (ja, jb, ja1, ipca, jb1)
      implicit none (type, external)
      integer, intent(in) :: ja
      integer, intent(in) :: jb
      integer, intent(in) :: ja1
      integer, intent(in) :: ipca
      integer, intent(in) :: jb1
      end subroutine
      end interface
      end module
