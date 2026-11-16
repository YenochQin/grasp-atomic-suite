      module couple_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, parameter :: mangm = 60
      integer, parameter :: mtriad = 12
      integer, dimension(mangm) :: j1
      integer, dimension(mtriad,3) :: j2, j3
      integer :: mja, nja
      logical, dimension(mangm) :: free
      end module couple_c
