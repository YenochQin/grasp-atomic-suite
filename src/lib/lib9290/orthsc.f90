module orthsc_i
    implicit none (type, external)
    private
    public :: orthsc
contains
!***********************************************************************
!                                                                      *
      subroutine orthsc
!                                                                      *
!   This routine Schmidt orthogonalises radial wavefunctions.          *
!                                                                      *
!   Call(s) to: [LIB92]: RINT.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 14 Oct 1992   *
!                                                                      *
! Normalization of the orbitals moved out of the inner loop
! XHH 1997.02.14
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only:  max_orbitals
      use debug_c
      use def_c
      use grid_c
      use orb_c
      use wave_c, only: orbital_origin_parameters, orbital_last_grid_index, radial_large_components, radial_small_components
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rint_i, only: rint
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j(max_orbitals)
      integer :: l, nakl, kount, mtp0, k, mtp, i
      real(real64) :: ovlap(max_orbitals)
      real(real64) :: eps, ovrlap, dnorm, factor
      logical :: changed
!-----------------------------------------------
!
!XHH
!
!   Set tabulated values of the radial wavefunction to zero
!   if they are less than EPS
!
      eps = 0.01_real64*accy
!
!   Determine the number of interesting overlaps
!
!XHH wasting time, no use at all
!      NOVL = 0
!      DO 2 K = 1,NW-1
!         NAKK = NAK(K)
!         DO 1 L = K+1,NW
!            IF (NAKK .EQ. NAK(L)) NOVL = NOVL+1
!    1    CONTINUE
!    2 CONTINUE
!
!      IF (NOVL .EQ. 0) RETURN
!
      do l = 2, orbital_count
!
         nakl = orbital_kappa(l)
         kount = 0
!XHH MTP0 introduced to count the maximum number of points during
!    orthogonalization of the L-th orbital to other orbitals
!    The logical variable changed is initialized to .F.

         mtp0 = orbital_last_grid_index(l)
         changed = .false.
!
         do k = 1, l - 1
!
            if (orbital_kappa(k) /= nakl) cycle
!XHH
            changed = .true.
!
!   Compute overlap
!
            ovrlap = rint(l,k,0)
!
!   Schmidt orthogonalise
!
            kount = kount + 1
            j(kount) = k
            ovlap(kount) = ovrlap
!
            orbital_origin_parameters(l) = orbital_origin_parameters(l) - ovrlap*orbital_origin_parameters(k)
            mtp = max(orbital_last_grid_index(l),orbital_last_grid_index(k))
            mtp0 = max(mtp0,orbital_last_grid_index(k))

            radial_large_components(:mtp,l) = radial_large_components(:mtp,l) - ovrlap*radial_large_components(:mtp,k)
            radial_small_components(:mtp,l) = radial_small_components(:mtp,l) - ovrlap*radial_small_components(:mtp,k)
         end do
!
!   Normalise
!
!XHH Use MTP0 to replace MTP and only when the orbital is changed.
!    This is in accordance with the original version which had the
!    normalization etc within the inner K loop.

         if (changed) then
            mtp = mtp0

            orbital_last_grid_index(l) = mtp
            dnorm = rint(l,l,0)
            factor = 1.0_real64/sqrt(dnorm)
!
            orbital_origin_parameters(l) = factor*orbital_origin_parameters(l)
            radial_large_components(2:mtp,l) = factor*radial_large_components(2:mtp,l)
            radial_small_components(2:mtp,l) = factor*radial_small_components(2:mtp,l)
!
!   Find new MF(L)
!
            mtp = mtp + 1
    7       continue
            mtp = mtp - 1
!cjb       print *, '    MTP    = ',    MTP
!          print *, '        L  = ',        L
!          print *, ' PF(MTP,L) = ', PF(MTP,L)
!cjb Subscript out of range for array pf
           if ( mtp >= 1 ) then
!cjb
            if (abs(radial_large_components(mtp,l)) < eps) then
               radial_large_components(mtp,l) = 0.0_real64
               radial_small_components(mtp,l) = 0.0_real64
               goto 7
            else
               orbital_last_grid_index(l) = mtp
            end if
!cjb endif
           end if
         end if
         if (.not.(ldbpr(3) .and. kount>0)) cycle

!XHH Moved ahead
!            ENDIF
!    8    CONTINUE
!
!   Print overlap information
!
!---------------------------------------------------------------
! debug 97.02.14 done
! Check the orthonormalisation
! KOUNT and OVLAP are destroyed
! should be commented out after seeing the values
!
!      KOUNT = 0
!      DO k = 1, L
!            IF (NAK(K) .EQ. NAKL) THEN
!               KOUNT = KOUNT+1
!               OVLAP(KOUNT) = RINT (L,K,0)
!               J(KOUNT) = K
!            ENDIF
!      ENDDO
!      WRITE (*,301)
!     :      (OVLAP(I),NP(L),NH(L),NP(J(I)),NH(J(I)),I = 1,KOUNT)
!---------------------------------------------------------------
         write (99, 301) (ovlap(i),orbital_principal_n(l),orbital_symmetry_labels(l),orbital_principal_n(j(i)), &
            & orbital_symmetry_labels(j(i)),i=1,kount)
!
      end do
!
      return
!
  301 format(1p,5(2x,1d10.3," = <",1i2,1a2,"|",1i2,1a2,">"))
      return
!
      end subroutine orthsc
end module orthsc_i
