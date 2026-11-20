!***********************************************************************
!                                                                      *
      integer function iqr (isubsh, icsf)
!                                                                      *
!   IQR is the occupation of subshell ISUBSH in CSF  ICSF.             *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use orbital_csf_basis,         only: csf_subshell_occupations
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: isubsh
      integer, intent(in) :: icsf
!-----------------------------------------------
!
!
            iqr = csf_subshell_occupations(isubsh,icsf)
!
      return
      end function iqr
