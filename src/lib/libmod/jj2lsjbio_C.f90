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
      use vast_kind_param, only:  double
!
      real(double), dimension(:), pointer :: rlev_ener_1,rlev_ener_2
      character(len=64), pointer, dimension(:) :: string_csf1
      character(len=64), pointer, dimension(:) :: string_csf2
      integer :: nvectoti,nvectotf,iopen_status1,iopen_status2
      integer, dimension(:), pointer :: lev_pos_1,lev_j_1,lev_par_1
      integer, dimension(:), pointer :: lev_pos_2,lev_j_2,lev_par_2
!
      end module jj2lsjbio_c
