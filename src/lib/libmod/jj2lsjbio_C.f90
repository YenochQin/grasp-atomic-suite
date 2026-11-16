!
!***********************************************************************
!                                                                      *
      module jj2lsjbio_c
!                                                                      *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2017   *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!
      real(real64), pointer :: rlev_ener_1(:),rlev_ener_2(:)
      character(len=64), pointer :: string_csf1(:)
      character(len=64), pointer :: string_csf2(:)
      integer :: nvectoti,nvectotf,iopen_status1,iopen_status2
      integer, pointer :: lev_pos_1(:),lev_j_1(:),lev_par_1(:)
      integer, pointer :: lev_pos_2(:),lev_j_2(:),lev_par_2(:)
!
      end module jj2lsjbio_c
