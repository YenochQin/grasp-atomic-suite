!***********************************************************************
!                                                                      *
      subroutine orthor(j, inv)
!
!                                                                      *
!   This routine Schmidt orthogonalizes orbital  J  to  all orbitals   *
!   which  have  better  self-consistency.  Note that fixed orbitals   *
!   have the best self-consistency.                                    *
!                                                                      *
!   Call(s) to: [LIB92]: COUNT, RINT.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!                                                                      *
!***********************************************************************
! INV not initialized here but used as
! INV = INV + 1
! A bug regarding IORDER fixed and new schemes introduced. See comments
! Anyway this routine was not used anywhere.
!XHH 1997.02.21
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:54:29   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_radial_points
      use def_c
      use grid_c
      use ovl_c
      use orb_c
      use orba_c
      use invt_c
      use scf_c
      use iounit_c
      use fixd_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rint_i
!      USE count_I
      use count_i, only: count
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      integer , intent(inout) :: inv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kraw, k, mtp, i, nncff, jgg
      real(real64) :: eps, ovrlap, dnorm, factor, sgn
      logical :: check, changed
!-----------------------------------------------
!
!XHH Added common /orba/ and /iounit/

!
!ww      EPS = ACCY*0.1D 00
      eps = accy*0.01_real64
!
      check = .not.noinvt(j)
!
!XHH
! Bug fixed. Here J is actually IORDER(J_raw), thus K should be
! treated the same way.
      changed = .false.
!      DO 4 K = 1,NW
      do kraw = 1, orbital_count
         k = iorder(kraw)
!         write(istde,*) '***',kraw,k,np(k),nh(k),scnsty(k),'***'
!
!XHH orbitals with higher self-consistency are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) .AND.
!     :        (SCNSTY(K) .LT. SCNSTY(J)) ) THEN
!XHH All orbitals are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) ) THEN
!XHH orbitals before the current and unchanged ones are considered
         if (.not.(orbital_kappa(k)==orbital_kappa(j) .and. (k<j .or. is_orbital_fixed(k)))) cycle

         changed = .true.
!
!   Compute overlap
!
         ovrlap = rint(j,k,0)
!
!   Schmidt orthogonalise
!
         orbital_origin_parameters(j) = orbital_origin_parameters(j) - ovrlap*orbital_origin_parameters(k)
         mtp = max(orbital_last_grid_index(j),orbital_last_grid_index(k))
         radial_large_components(:mtp,j) = radial_large_components(:mtp,j) - ovrlap*radial_large_components(:mtp,k)
         radial_small_components(:mtp,j) = radial_small_components(:mtp,j) - ovrlap*radial_small_components(:mtp,k)
      end do
!XHH
      if (changed) then
!
!   Normalise
!
         orbital_last_grid_index(j) = mtp
         dnorm = rint(j,j,0)
         factor = 1.0_real64/sqrt(dnorm)
!
!   Determine if inversion is necessary
!
         if (check) then
            call count (radial_large_components(:max_radial_points,j), mtp, nncff, sgn)
            if (sgn < 0.0_real64) then
               inv = inv + 1
               factor = -factor
            endif
         endif
!
!   Perform normalization and/or inversion
!
         orbital_origin_parameters(j) = factor*orbital_origin_parameters(j)
         radial_large_components(2:mtp,j) = factor*radial_large_components(2:mtp,j)
         radial_small_components(2:mtp,j) = factor*radial_small_components(2:mtp,j)
!
!   Find new MF(J)
!
         mtp = mtp + 1
    3    continue
         mtp = mtp - 1
         if (abs(radial_large_components(mtp,j)) < eps) then
            radial_large_components(mtp,j) = 0.0_real64
            radial_small_components(mtp,j) = 0.0_real64
            go to 3
         else
            orbital_last_grid_index(j) = mtp
         endif
!
!XHH two lines moved forward and added an ENDIF
!         ENDIF
!    4 CONTINUE
      endif
!
      return
      end subroutine orthor
