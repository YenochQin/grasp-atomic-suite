!
!***********************************************************************
!                                                                      *
      module syma_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double
      integer, parameter::nblk0 = 50    ! Maximum number of blocks
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, dimension(:), pointer :: iatjpo, iaspar
      integer, dimension(:), pointer :: iatjpoff, iasparff
      integer, dimension(:), pointer :: iatjpoii, iasparii
      integer, dimension(nblk0) :: jpgg
      end module syma_c
