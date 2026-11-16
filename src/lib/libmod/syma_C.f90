!
!***********************************************************************
!                                                                      *
      module syma_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer, parameter::nblk0 = 50    ! Maximum number of blocks
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, pointer :: iatjpo(:), iaspar(:)
      integer, pointer :: iatjpoff(:), iasparff(:)
      integer, pointer :: iatjpoii(:), iasparii(:)
      integer :: jpgg(nblk0)
      end module syma_c
