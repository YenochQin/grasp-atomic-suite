!
!***********************************************************************
!                                                                      *
      module def_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double
      use parameter_def,   only:  nnnp
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double) :: tenmax, expmax, expmin, precis
      real(double) :: aucm, auev, ccms, fasi, fbsi
      real(double) :: fmtoau, aumamu, b1
      integer :: ioncty, nelec
      real(double) :: emn, z
      integer :: ionctyff, nelecff
      real(double) :: emnff, zff
      integer :: ionctyii, nelecii
      real(double) :: emnii, zii
      integer :: nelecr
      real(double) :: c
      real(double) :: empam, rbcm
      integer :: nscf, nsic, nsolv
      real(double) :: accy
      real(double), dimension(:), pointer :: wt, weight
      integer :: ncmin, ncmax
      integer, dimension(:), pointer :: iccmin
      real(double) :: cvac, pi
      real(double), dimension(nnnp) :: dp, dq
      end module def_c
