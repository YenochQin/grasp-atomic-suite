!***********************************************************************
!                                                                      *
      subroutine bessel(ia, ib, ik, iw, k)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!                                                                      *
!   This routine evaluates the functions                               *
!                                                                      *
!                  (2K+1)!!                                            *
!      BESSJ  =    --------  J  (w   r) - 1  =  PHI  (w   r) - 1       *
!                         K   K   ab               K   ab              *
!                  (w   r)                                             *
!                    ab                                                *
!   and                                                                *
!                         K+1                                          *
!                  (w   r)                                             *
!                    ab                                                *
!      BESSN  =  - --------  N  (w   r) - 1  =  PSI  (w   r) - 1       *
!                  (2K-1)!!   K   ab               K   ab              *
!                                                                      *
!                                                                      *
!   where J and N ARE spherical Bessel functions, and                  *
!                                                                      *
!                w   = ABS ( E(IA) - E(IB) )/c                         *
!                 ab                                                   *
!                                                                      *
!   where  E(I)  is the eigenvalue for orbital  I . The writeup (B J   *
!   McKenzie, I P Grant, and P H Norrington, Computer Phys Commun 21   *
!   (1980) 233-246) is incorrect in its description of the output of   *
!   this routine.                                                      *
!                                                                      *
!   The  routine uses equations given in M Abramowitz and I A STegun   *
!   to evaluate the functions. Devices are used to reduce the number   *
!   of actual evaluations of these functions.                          *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M O D U L E S
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use bess_c
      use debug_c
      use def_c
      use grid_c
      use orb_c
      use stor_c
      use wfac_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      integer, intent(in) :: ik
      integer, intent(in) :: iw
      integer, intent(in) :: k
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icode, iwkp, ikkp, l, i, nn, j, jchan, irem, iswap
      real(real64) :: epsi, w, wa, xbess1, xbess2, s1, s2, dfnm, dfn, ssn, scn&
         , sn, cn, obwa, b, skeep
!-----------------------------------------------
!
      epsi = dsqrt(0.1_real64*accy)
!
!   Form unique label symmetric in IA, IB
!
      icode = max(ia,ib) + key*(min(ia,ib) + key*k)
!
!   Function in position; return
!
      if (icode == keep(ik,iw)) return
!
!   Function not in position; is it available in BESS arrays?
!
      w = wfact*dabs(orbital_energy_parameters(ia)-orbital_energy_parameters(ib))/speed_of_light_au
      wij(iw) = w
!
      do iwkp = 1, 2
         do ikkp = 1, 2
            if (keep(ikkp,iwkp) /= icode) cycle
!
!   Function found move into position
!
            keep(ik,iw) = icode
            if (ldbpr(7)) write (99, 302) orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib) &
               & , orbital_symmetry_labels(ib), k, &
               ikkp, iwkp, ik, iw
            bessj(ik,iw,:radial_point_count) = bessj(ikkp,iwkp,:radial_point_count)
            bessn(ik,iw,:radial_point_count) = bessn(ikkp,iwkp,:radial_point_count)
            return
         end do
      end do
!
!   Function not found; evaluate it
!
      if (ldbpr(7)) write (99, 303) orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
         & orbital_symmetry_labels(ib), k, ik, iw
!
      keep(ik,iw) = icode
!
      if (w < epsi**2) then
!
!   Negligible w
!
         bessj(ik,iw,:radial_point_count) = 0.0_real64
         bessn(ik,iw,:radial_point_count) = 0.0_real64
         return
!
      endif
!
      nn = k
!
      bessj(ik,iw,1) = 0.0_real64
      bessn(ik,iw,1) = 0.0_real64
!
!   Use a four-term power series for low w*r
!
      l5: do j = 2, radial_point_count
         wa = -0.5_real64*(radial_grid(j)*w)**2
         xbess1 = 1.0_real64
         xbess2 = 1.0_real64
         s1 = 0.0_real64
         s2 = 0.0_real64
         do i = 1, 4
            xbess1 = xbess1*wa/dble(i*(2*(nn + i) + 1))
            xbess2 = xbess2*wa/dble(i*(2*(i - nn) - 1))
            s1 = s1 + xbess1
            s2 = s2 + xbess2
            if (dabs(xbess1)>=dabs(s1)*epsi .or. dabs(xbess2)>=dabs(s2)*epsi) &
               cycle
            bessj(ik,iw,j) = s1
            bessn(ik,iw,j) = s2
            cycle  l5
         end do
         jchan = j
         go to 6
      end do l5
!
!   If here then calculated whole array using four-term power
!   series.  Hence return
!
      return
!
!   Use sin/cos expansion when power series requires more than
!   four terms terms to converge
!
    6 continue
      if (nn == 0) then
         dfnm = 1.0_real64
         dfn = 1.0_real64
      else
         dfnm = 1.0_real64
         do i = 3, 2*nn - 1, 2
            dfnm = dfnm*dble(i)
         end do
         dfn = dfnm*dble(2*nn + 1)
      endif
      dfnm = 1.0_real64/dfnm
!
      irem = mod(nn,4)
!
      select case (irem)
      case (1)
!
!   NN = 1, 5, 9, ...
!
         ssn = -1.0_real64
         scn = 1.0_real64
         iswap = 1
!
      case (2)
!
!   N = 2, 6, 10, ....
!
         ssn = -1.0_real64
         scn = -1.0_real64
         iswap = 0
!
      case (3)
!
!   NN = 3, 7, 11,...
!
         ssn = 1.0_real64
         scn = -1.0_real64
         iswap = 1
!
      case default
!
!   NN = 0, 4, 8,...
!
         ssn = 1.0_real64
         scn = 1.0_real64
         iswap = 0
!
      end select
!
      do j = jchan, radial_point_count
         wa = w*radial_grid(j)
         if (iswap == 0) then
            sn = ssn*dsin(wa)
            cn = scn*dcos(wa)
         else
            sn = ssn*dcos(wa)
            cn = scn*dsin(wa)
         endif
         obwa = 1.0_real64/wa
         b = obwa
         s1 = b*sn
         s2 = b*cn
         do i = 1, nn
            skeep = sn
            sn = cn
            cn = -skeep
            b = b*obwa*dble((nn + i)*(nn - i + 1))/dble(2*i)
            s1 = s1 + b*sn
            s2 = s2 + b*cn
         end do
         s1 = s1*dfn/wa**nn - 1.0_real64
         s2 = s2*wa**(nn + 1)*dfnm - 1.0_real64
         bessj(ik,iw,j) = s1
         bessn(ik,iw,j) = s2
      end do
      return
!
  303 format(93x,i2,a2,2x,i2,a2,2x,i2,2x,'New',6x,'(',i2,',',i2,')')
  302 format(93x,i2,a2,2x,i2,a2,2x,i2,2x,'(',i2,',',i2,')',2x,'(',i2,',',i2,')'&
         )
      return
!
      end subroutine bessel
