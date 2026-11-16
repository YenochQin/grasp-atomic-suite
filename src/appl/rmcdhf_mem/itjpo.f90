!***********************************************************************
!                                                                      *
      integer function itjpo (icsf)
!                                                                      *
!   ITJPO is the value of 2J+1 for CSF number ICSF.                    *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:45   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use stat_c,        only: csf_coupling_quantum_codes
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: icsf
!-----------------------------------------------
      itjpo = csf_coupling_quantum_codes(max_orbitals,icsf)
      if (itjpo > 127) itjpo = 256 - itjpo
      itjpo = iabs (itjpo)
      return
      end function itjpo
