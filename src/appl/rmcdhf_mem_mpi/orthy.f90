!***********************************************************************
!
      subroutine orthy(nw, jp, lsort)
!
! nw     - input, total number of orbitals
! jp     - input, orbital jp, usually updated from  solve. It is
!          the one after IORDER
! lsort  - input, logical, .t. then sort according in decreasing
!          self-consistency order.
!
! The following quantities are not used outside this routine
!
! kfixed - output, the number of fixed orbitals of the same kappa
!          (specified by nak(jp)) as the orbital jp.
! ktotal - output, the number of all orbitals of the same kappa as the
!          orbital jp.
! kindx  - output, an index array containing positions of the orbitals
!          in the decreasing self-consistency order. It has the size of
!          ktotal.
!
! Orbitals are grouped as fixed, unfixed (non-correlation, correlation).
! .For i=1 to kfixed, kindx(i) has the same order as array iorder (though
!  this can be modified by requiring the smallest principal quantum number
!  come first).
! .For i=kfixed+1 to knon, kindx(i) is ordered as SCN/E**2.
! .For i=knon+1 to ktotal, kindx(i) is ordered as SCN.
! .When there is degeneracy (same SCN/E**2 or same SCN), iorder is used.
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:57:18   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_orbitals
      use def_c
      use orb_c, only: orbital_kappa, orbital_energy_parameters
      use orba_c, only: iorder
      use corre_c, only: lcorre
      use fixd_c, only: is_orbital_fixed
      use invt_c, only: noinvt
      use iounit_c
      use scf_c, only: scnsty
      use wave_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rint_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nw
      integer , intent(in) :: jp
      logical , intent(in) :: lsort
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j1
      integer  :: kindx(max_orbitals)
      integer :: kfixed, nakj, ktotal, i, ki, j, knon, ishift, kj, lraw, l, &
         nakl, mtp0, kraw, k, mtp
      real(real64) :: eps, ovrlap, dnorm, factor
      logical :: check
!-----------------------------------------------
!
! For inverse - taken from orthor.f which was not used anywhere
! Not used currently
!-----------------------------------------------------------------------
      check = .not.noinvt(jp)

      eps = 0.01_real64*accy
      kfixed = 0
      nakj = orbital_kappa(jp)

!-----------------------------------------------------------------------
! Find orbitals of the same kappa as orbital jp.
!-----------------------------------------------------------------------
      ktotal = 0
      do i = 1, nw
!         Orbital jp is supposed to be an ordered one, so should be k
         ki = iorder(i)
         if (orbital_kappa(ki) /= nakj) cycle
         ktotal = ktotal + 1
         kindx(ktotal) = ki
      end do

!-----------------------------------------------------------------------
! Find fixed orbitals and separate them from the rest.
!-----------------------------------------------------------------------

!      Place the fixed orbitals at the begining of the array, no
!      sorting is done here. And find the number of non-correlation
!      orbitals.

      j = 0
      knon = 0
      do i = 1, ktotal
         ki = kindx(i)
         if (is_orbital_fixed(ki)) then

            j = j + 1

             ! Shift array elements rather than exchanging
            kindx(i:1+j:(-1)) = kindx(i-1:j:(-1))
!            kindx(i) = kindx(j)

            kindx(j) = ki

         else if (.not.lcorre(ki)) then
            knon = knon + 1
         endif
      end do
      kfixed = j

      if (lsort) then

!-----------------------------------------------------------------------
! Separate correlation from non-correlation orbitals
!-----------------------------------------------------------------------

         j = kfixed
         do i = kfixed + 1, ktotal
            ki = kindx(i)
            if (lcorre(ki)) cycle

            j = j + 1

                ! Shift array elements rather than exchanging
            kindx(i:1+j:(-1)) = kindx(i-1:j:(-1))
!               kindx(i) = kindx(j)

            kindx(j) = ki
         end do

!-----------------------------------------------------------------------
! Sort non-correlation and correlation orbitals separately
! This can be done indepently since knon is known
!-----------------------------------------------------------------------

!           non-correlation, using criteria scnsty/E^2

         do i = kfixed + 1, kfixed + knon
            ki = kindx(i)
            do j = i + 1, kfixed + knon
               kj = kindx(j)
               if (scnsty(kj)/orbital_energy_parameters(kj)**2 >= scnsty(ki)/orbital_energy_parameters(ki)**2) cycle

                   ! No need to do shifting, exchanging is fine
               kindx(j) = ki
               ki = kj
            end do
            kindx(i) = ki
         end do

!           correlation orbitals, using criteria scnsty

         do i = kfixed + knon + 1, ktotal
            ki = kindx(i)
            do j = i + 1, ktotal
               kj = kindx(j)
               if (scnsty(kj) >= scnsty(ki)) cycle

                   ! No need to do shifting, exchanging is fine
               kindx(j) = ki
               ki = kj
            end do
            kindx(i) = ki
         end do

      endif

! Finished sorting.
! Schmidt orthogonalize all orbitals of the same kappa
! The fixed orbitals are not changed

      do lraw = kfixed + 1, ktotal
         l = kindx(lraw)

         nakl = orbital_kappa(l)

         mtp0 = orbital_last_grid_index(l)

         do kraw = 1, lraw - 1
            k = kindx(kraw)
            ovrlap = rint(l,k,0)

!   Schmidt orthogonalise

            orbital_origin_parameters(l) = orbital_origin_parameters(l) - ovrlap*orbital_origin_parameters(k)
            mtp = max(orbital_last_grid_index(l),orbital_last_grid_index(k))
            mtp0 = max(mtp0,orbital_last_grid_index(k))

            radial_large_components(:mtp,l) = radial_large_components(:mtp,l) - ovrlap*radial_large_components(:mtp,k)
            radial_small_components(:mtp,l) = radial_small_components(:mtp,l) - ovrlap*radial_small_components(:mtp,k)
         end do

!   Normalise

         mtp = mtp0

         orbital_last_grid_index(l) = mtp
         dnorm = rint(l,l,0)
         factor = 1._real64/sqrt(dnorm)

         if (orbital_origin_parameters(l) < 0._real64) factor = -factor

         orbital_origin_parameters(l) = factor*orbital_origin_parameters(l)
         radial_large_components(2:mtp,l) = factor*radial_large_components(2:mtp,l)
         radial_small_components(2:mtp,l) = factor*radial_small_components(2:mtp,l)

!   Find new MF(L)

         mtp = mtp + 1
   20    continue
         mtp = mtp - 1
         if (abs(radial_large_components(mtp,l)) < eps) then
            radial_large_components(mtp,l) = 0._real64
            radial_small_components(mtp,l) = 0._real64
            go to 20
         else
            orbital_last_grid_index(l) = mtp
         endif

      end do

      return
      end subroutine orthy
