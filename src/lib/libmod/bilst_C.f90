      module bilst_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(6) :: ndtpa, ntpi
      integer, dimension(:), pointer :: indtp1, indtp2, indtp3, indtp4, &
                                        indtp5, indtp6
      real(double), dimension(:), pointer :: valtp1, valtp2, valtp3, valtp4, &
                                        valtp5, valtp6
      logical, dimension(6) :: first
      end module bilst_c
