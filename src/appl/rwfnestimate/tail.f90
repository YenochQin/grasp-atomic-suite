

!***********************************************************************
!                                                                      *
      subroutine tail(iorb, p, q, jp, mtp)
!                                                                      *
!   This subroutine begins the inward integration of the homogeneous   *
!   Dirac radial equation. With only minor modifications, the series   *
!   given by  J E Sienkiewicz  and W E Baylis, J Phys B: At Mol Phys   *
!   20 (1987) 5145-5156, p 5155, is used.                              *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def, only: max_radial_points
      use def_c, only: accy, speed_of_light_au, expmax
      use grid_c
      use orb_c
      use pote_c, only: yp
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iorb
      integer , intent(in) :: jp
      integer , intent(out) :: mtp
      real(real64) , intent(inout) :: p(max_radial_points)
      real(real64) , intent(out) :: q(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nm4, m, loc
      real(real64) :: eps, bige, bigebc, bgebc2, beta, fk, rjp, qem, offset, t&
         , zloc, gam2, zlbb, dnu, fkmzbb, sump, sumq, em, cm, emfact, ovltrm, &
         pterm, qterm, exptrm
!-----------------------------------------------
!
!
!   Initialize
!
      eps = accy*0.1_real64
      bige = -orbital_energy_parameters(iorb)
      bigebc = bige/speed_of_light_au
      bgebc2 = bigebc/speed_of_light_au
      beta = sqrt((-bige*(2.0_real64 + bgebc2)))
      fk = dble(orbital_kappa(iorb))
!
!   Find MTP
!
      i = jp
      nm4 = radial_point_count - 4
      rjp = radial_grid(jp)
      qem = 0.25_real64*expmax
    1 continue
      i = i + 1
      if (i <= nm4) then
         if (beta*(radial_grid(i)-rjp) > qem) then
            mtp = i
         else
            go to 1
         endif
      else
         write (*, 300)
         mtp = nm4
      endif
!
!   Compute offset for exponential function
!
      offset = beta*radial_grid(mtp)
!
!   Tabulate tail points
!
      do i = mtp, radial_point_count
!
         t = 2.0_real64*beta*radial_grid(i)
!
         zloc = yp(i)
         gam2 = fk**2 - (zloc/speed_of_light_au)**2
         zlbb = zloc/beta
         dnu = zlbb*(1.0_real64 + bgebc2)
         fkmzbb = fk - zlbb
!
         m = -1
         sump = 0.0_real64
         sumq = 0.0_real64
         m = m + 1
         em = dble(m)
         if (m == 0) then
            cm = 1.0_real64
            emfact = 1.0_real64
         else
            cm = cm*(gam2 - (dnu - em)**2)
            emfact = emfact*em
         endif
         ovltrm = cm*(t**(dnu - em)/emfact)
         pterm = ovltrm*(fkmzbb + em)*beta
         qterm = ovltrm*(fkmzbb - em)*bigebc
         sump = sump + pterm
         sumq = sumq + qterm
         do while(abs(pterm/sump)>=eps .or. abs(qterm/sumq)>=eps)
            m = m + 1
            em = dble(m)
            if (m == 0) then
               cm = 1.0_real64
               emfact = 1.0_real64
            else
               cm = cm*(gam2 - (dnu - em)**2)
               emfact = emfact*em
            endif
            ovltrm = cm*(t**(dnu - em)/emfact)
            pterm = ovltrm*(fkmzbb + em)*beta
            qterm = ovltrm*(fkmzbb - em)*bigebc
            sump = sump + pterm
            sumq = sumq + qterm
         end do
         exptrm = exp((-0.5_real64*t) + offset)
         p(i) = sump*exptrm
         q(i) = sumq*exptrm
         if (p(i) /= 0.0_real64) cycle
         loc = i + 1
         go to 4
      end do
      loc = radial_point_count + 1
!
    4 continue
      p(loc:radial_point_count) = 0.0_real64
      q(loc:radial_point_count) = 0.0_real64
!
      return
!
  300 format('TAIL: Grid may be of insufficient extent')
      return
!
      end subroutine tail
