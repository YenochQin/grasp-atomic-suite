!
!***********************************************************************
!                                                                      *
      module stat_c
!
!***********************************************************************
!                                                                      *
      use vast_kind_param, only:  double, byte
!...Created by Pacific-Sierra Research 77to90  4.3E  07:21:55   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(:,:,:), pointer :: jqsar
      integer(byte), dimension(:,:,:), pointer :: jqsa
      integer, dimension(:,:), pointer :: jcupar
      integer(byte), dimension(:,:), pointer :: jcupa
      end module stat_c
