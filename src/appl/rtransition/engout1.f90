!***********************************************************************
!                                                                      *
      subroutine engout1(eav, e, jtot, ipar, ilev, nn, mode, k)
!                                                                      *
!   This  subroutine prints  energy levels, splittings, and energies   *
!   relative to the lowest in  Hartrees, Kaysers, and  eV, using the   *
!   reduced mass corrected value for the Rydberg. If MODE is 0, only   *
!   the eigenenergies are printed. If  MODE  is 1, the eigenenergies   *
!   and separations are printed. If MODE is 2, the eigenenergies and   *
!   energies relative to level 1 are printed. If MODE is 3, the eig-   *
!   enenergies,  separations,  and energies relative to level  1 are   *
!   printed.                                                           *
!                                          Last updated: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c, only: aucm, auev, ccms, fasi, fbsi
      use jlabl_c, labj=>jlbr, labp=>jlbp
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nn
      integer  :: mode
      integer , intent(in) :: k
      real(double) , intent(in) :: eav
      integer , intent(in) :: jtot(nn)
      integer , intent(in) :: ipar(nn)
      integer , intent(in) :: ilev(nn)
      real(double) , intent(in) :: e(nn)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, ip
      real(double) :: eau, ecm, eev
!-----------------------------------------------
!
!   Always print the eigenenergies
!
      if (k == 1) write (24, 299)
      if (k == 2) write (24, 300)
      write (24, 301)
      do j = 1, nn
         i = ilev(j)
         eau = e(j) + eav
         ecm = eau*aucm
         eev = eau*auev
         ip = (ipar(j)+3)/2
         write (24, 302) i, labj(jtot(j)), labp(ip), eau, ecm
      end do
!
      return
!
  299 format('Eigenenergies for the initial state list')
  300 format('Eigenenergies for the final state list')
  301 format('Level  J Parity',10x,'Hartrees',18x,'Kaysers')
  302 format(1i3,4x,2a4,1p,2d25.15)
  303 format('Energy of each level relative to immediately lower',' level:')
  304 format('Energy of each level relative to lowest level:')
      return
!
      end subroutine engout1
