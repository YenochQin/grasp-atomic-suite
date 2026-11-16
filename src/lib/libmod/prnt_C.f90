!
!***********************************************************************
!                                                                      *
      module prnt_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double, byte
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: nvec, nvecmx
      integer :: nvecff, nvecmxff
      integer :: nvecii, nvecmxii
      real(double) :: pnivecii
      integer, dimension(:), pointer :: ivec, ivecff, ivecii
      end module prnt_c
