!
!***********************************************************************
!                                                                      *
      module orb_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double, byte
      use parameter_def,   only:  nnnw
!...Created by Pacific-Sierra Research 77to90  4.3E  08:57:22  12/25/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17

      character(len=2), dimension(nnnw) :: nh
      character(len=2), dimension(nnnw) :: nhr
      real(double), dimension(nnnw) :: e, gama, ped
      integer :: ncf, nw, ncfr, nwr
      integer(byte), dimension(:,:), pointer :: iqa
      integer, dimension(:,:), pointer :: iqar
      integer, dimension(nnnw) :: np, nak, npr, nakr
      integer, dimension(nnnw) :: nkl, nkj
      end module orb_c
