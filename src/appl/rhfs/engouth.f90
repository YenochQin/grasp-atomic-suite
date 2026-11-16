!***********************************************************************
!                                                                      *
      subroutine engouth(eav, e, jtot, ipar, ilev, nn, mode)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c, only: hartree_to_wavenumber_cm, hartree_to_ev, speed_of_light_cm_per_s, einstein_a_au_to_si, einstein_b_au_to_si
      use jlabl_c, labj=> jlbr, labp=>jlbp
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nn
      integer, intent(in) :: mode
      real(real64), intent(in) :: eav
      integer, intent(in) :: jtot(nn)
      integer, intent(in) :: ipar(nn)
      integer, intent(in) :: ilev(nn)
      real(real64), intent(in) :: e(nn)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: j, i, ip
      real(real64) :: eau, ecm, eev
!-----------------------------------------------
!
!   Always print the eigenenergies
!
      write (24, 300)
      write (24, 301)
      write (29, 300)
      write (29, 301)
      do j = 1, nn
         i = ilev(j)
         eau = e(j) + eav
         ecm = eau*hartree_to_wavenumber_cm
         eev = eau*hartree_to_ev
         ip = (ipar(j)+3)/2
         write (24, 302) i, labj(jtot(j)), labp(ip), eau, ecm, eev
         write (29, 302) i, labj(jtot(j)), labp(ip), eau, ecm, eev
      end do
!
      if (nn > 1) then
!
!   Energy separations
!
         if (mode==1 .or. mode==3) then
            write (24, 303)
            write (24, 301)
            write (29, 303)
            write (29, 301)
            do j = 2, nn
               i = ilev(j)
               eau = e(j) - e(j-1)
               ecm = eau*hartree_to_wavenumber_cm
               eev = eau*hartree_to_ev
               ip = (ipar(j)+3)/2
               write (24, 302) i, labj(jtot(j)), labp(ip), eau, ecm, eev
               write (29, 302) i, labj(jtot(j)), labp(ip), eau, ecm, eev
            end do
         endif
!
!   Energies relative to level 1
!
         if (mode==2 .or. mode==3) then
            write (24, 304)
            write (24, 301)
            write (29, 304)
            write (29, 301)
            do j = 2, nn
               i = ilev(j)
               eau = e(j) - e(1)
               ecm = eau*hartree_to_wavenumber_cm
               eev = eau*hartree_to_ev
               ip = (ipar(j)+3)/2
               write (24, 302) i, labj(jtot(j)), labp(ip), eau, ecm, eev
               write (29, 302) i, labj(jtot(j)), labp(ip), eau, ecm, eev
            end do
         endif
!
      endif
!
  300 format(/,'Eigenenergies:')
  301 format(/,'Level  J Parity',7x,'Hartrees',14x,'Kaysers',16x,'eV'/)
  302 format(1i3,2x,2a4,1p,3d22.14)
  303 format(/,'Energy of each level relative to immediately lower',' level:')
  304 format(/,'Energy of each level relative to lowest level:')
      return
!
      end subroutine engouth
