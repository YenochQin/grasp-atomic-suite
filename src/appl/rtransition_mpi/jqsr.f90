!***********************************************************************
!                                                                      *
      integer function jqsr (iwhich, isubsh, icsf)
!                                                                      *
!   JQSR is a subshell quantum number for subshell ISUBSH in configu-  *
!   ration state function  ICSF:  the seniority if IWHICH is 1;  the   *
!   quantum number w if IWHICH is 2, and 2J+1 if IWHICH is 3.          *
!                                                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use stat_c,        only: jqsar
      use orb_c,         only: ncfr, nwr
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: iwhich
      integer, intent(in) :: isubsh
      integer, intent(in) :: icsf
!-----------------------------------------------
!
             jqsr = jqsar(isubsh,iwhich,icsf)
!
      return
      end function jqsr
