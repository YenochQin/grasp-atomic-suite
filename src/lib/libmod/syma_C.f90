!
!***********************************************************************
!                                                                      *
      module syma_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer, parameter::nblk0 = 50    ! Maximum number of blocks
      integer, pointer :: iatjpo(:), iaspar(:)
      integer, pointer :: iatjpoff(:), iasparff(:)
      integer, pointer :: iatjpoii(:), iasparii(:)
      integer :: jpgg(nblk0)
      end module syma_c
