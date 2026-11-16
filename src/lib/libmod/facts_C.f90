!
!***********************************************************************
!                                                                      *
      module facts_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  10:46:52   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, parameter :: mfact = 500
      real(real64) :: gam(mfact)
      end module facts_c
