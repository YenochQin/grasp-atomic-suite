!
!***********************************************************************
!                                                                      *
      module npar_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: nuclear_model_parameter_count
      real(real64) :: nuclear_model_parameters(2)
      end module npar_c

      module nsmdat_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:36:34  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(real64) :: sqn,dmomnm,qmomb
      real(real64) :: hfsi, hfsd, hfsq
      real(real64) :: smsi, smsd, smsq
      end module nsmdat_c
