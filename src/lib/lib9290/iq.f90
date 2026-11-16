!***********************************************************************
!                                                                      *
      integer function iq (isubsh, icsf)
!                                                                      *
!   IQ is the occupation of subshell ISUBSH in CSF  ICSF.              *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 30 Oct 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:38   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use orb_c,         only: csf_subshell_occupations, csf_count
      use iounit_c,      only: stderr_unit
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: isubsh
      integer, intent(in) :: icsf
!-----------------------------------------------
!
      iq = csf_subshell_occupations(isubsh,icsf)
!
      return
      end function iq
