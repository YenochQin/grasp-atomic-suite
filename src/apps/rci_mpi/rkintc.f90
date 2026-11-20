!***********************************************************************
!                                                                      *
      subroutine rkintc(orbital_a, orbital_b, orbital_c, orbital_d, multipole_rank, tegral)
!                                                                      *
!                         k                                            *
!   This routine returns R (abcd) integrals.                           *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use distributed_slater_integrals, only: indteirk, valteirk, rk_local_count
      use orbital_csf_basis, only: orbital_count
      use coulomb_rank_offsets, only: kmax, kstart
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: orbital_a, orbital_b, orbital_c, orbital_d
      integer, intent(in) :: multipole_rank
      real(real64), intent(out) :: tegral
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: key, iswap, index, jl, ju, jm, loc
      logical :: found, first
!-----------------------------------------------
!
      key = orbital_count + 1
!
!   Ensure that the indices are in `canonical' order
!   Compute the composite (packed) index
!


      if (orbital_a > orbital_c) then
         iswap = orbital_c
         orbital_c = orbital_a
         orbital_a = iswap
      end if
      if (orbital_b > orbital_d) then
         iswap = orbital_d
         orbital_d = orbital_b
         orbital_b = iswap
      end if
      if (orbital_a > orbital_b) then
         iswap = orbital_b
         orbital_b = orbital_a
         orbital_a = iswap
         iswap = orbital_d
         orbital_d = orbital_c
         orbital_c = iswap
      end if

      index = ((orbital_a*key + orbital_b)*key + orbital_c)*key + orbital_d
!
      jl = kstart(multipole_rank)
      ju = kstart(multipole_rank+1) - 1
      if (index<indteirk(jl) .or. index>indteirk(ju)) then
         write (*, *) "Something wrong in rkintc"
         stop
      end if
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
    1 continue
      if (ju - jl > 1) then
         jm = (ju + jl)/2
         if (indteirk(jm) > index) then
            ju = jm
         else
            jl = jm
         end if
         goto 1
      end if
!
!   The range is bracketed to the extent possible
!
      if (index == indteirk(ju)) then
         loc = ju
      else if (index == indteirk(jl)) then
         loc = jl
      else
         write (*, *) multipole_rank, orbital_a, orbital_b, orbital_c, orbital_d, index
         write (*, *) "Rkintc Integral not found"
         stop
      end if
!
!   Return the value of the integral
!   from storage

      tegral = valteirk(loc)
!
      return
      end subroutine rkintc
