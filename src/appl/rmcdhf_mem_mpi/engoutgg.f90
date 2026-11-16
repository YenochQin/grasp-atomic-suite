!***********************************************************************
!                                                                      *
      subroutine engoutgg(e,ilev,nn,mode)
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
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
      use jlabl_c, labj=>jlbr, labp=>jlbp
      use blkidx_c, only: double, idxblk
      use peav_c, only: double, eavblk
      use syma_c,          only: jpgg
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nn
      integer, intent(in) :: mode
      integer, intent(in) :: ilev(nn)
      real(real64), intent(in) :: e(nn)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jblock, i, ip, jtot
      real(real64) :: eav, eau, ecm, eev
!-----------------------------------------------
!
!
!   Always print the eigenenergies
!
      write (24, 300)
      write (24, 301)
      do j = 1, nn
         jblock = idxblk(j)
         eav = eavblk(jblock)
         i = ilev(j)
         eau = e(j) + eav
         ecm = eau*hartree_to_wavenumber_cm
         eev = eau*hartree_to_ev
!GG         IP = (IPAR(J)+3)/2
         jtot = iabs(jpgg(jblock))
         if(jpgg(jblock) >= 0) then
            ip = 2
         else
            ip = 1
         end if
         write (24, 302) i, labj(jtot), labp(ip), eau, ecm, eev
!GG         WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
      end do
!
      if (nn > 1) then
!
!   Energy separations
!
         if (mode==1 .or. mode==3) then
            write (24, 303)
            write (24, 301)
            do j = 2, nn
               i = ilev(j)
               eau = e(j) - e(j-1)
               ecm = eau*hartree_to_wavenumber_cm
               eev = eau*hartree_to_ev
!GG               IP = (IPAR(J)+3)/2
               jblock = idxblk(j)
               jtot = iabs(jpgg(jblock))
               if(jpgg(jblock) >= 0) then
                  ip = 2
               else
                  ip = 1
               end if
               write (24, 302) i, labj(jtot), labp(ip), eau, ecm, eev
!GG               WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
            end do
         end if
!
!   Energies relative to level 1
!
         if (mode==2 .or. mode==3) then
            write (24, 304)
            write (24, 301)
            do j = 2, nn
               i = ilev(j)
               eau = e(j) - e(1)
               ecm = eau*hartree_to_wavenumber_cm
               eev = eau*hartree_to_ev
!GG               IP = (IPAR(J)+3)/2
               jblock = idxblk(j)
               jtot = iabs(jpgg(jblock))
               if(jpgg(jblock) >= 0) then
                  ip = 2
               else
                  ip = 1
               end if
               write (24, 302) i, labj(jtot), labp(ip), eau, ecm, eev
!GG               WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
            end do
         end if
!
      end if
!
      return
!
  300 format(/,"Eigenenergies:")
  301 format(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)
  302 format(1i3,4x,2a4,1p,3d21.12)
  303 format(/,"Energy of each level relative to immediately lower"," level:")
  304 format(/,"Energy of each level relative to lowest level:")
      return
!
      end subroutine engoutgg
