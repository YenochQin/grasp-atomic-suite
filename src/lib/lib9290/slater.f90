module slater_i
    implicit none (type, external)
    private
    public :: slater
contains
!***********************************************************************
!                                                                      *
      real(real64) function slater (ia, ib, ic, id, k)
!                                                                      *
!   The value of this  function is the Slater integral                 *
!                                                                      *
!                               k                                      *
!                              R (abcd)                                *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD, YZK.                                    *
!                                                                      *
!                                         Last revision: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_c
      use grid_c
      use orb_c
      use tatb_c,          only: ta, tb, mtp
      use wave_c,          only: orbital_last_grid_index, radial_small_components, radial_large_components
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yzk_i, only: yzk
      use quad_i, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      integer, intent(in) :: ic
      integer, intent(in) :: id
      integer, intent(in) :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: result
!-----------------------------------------------
!
!
      call yzk (k, ib, id)
!
!   Multiply by second term, and obtain result by integration
!
      if (k==0 .and. ib==id) then
         mtp = min(orbital_last_grid_index(ia),orbital_last_grid_index(ic))
      else
         mtp = min(min(mtp,orbital_last_grid_index(ia)),orbital_last_grid_index(ic))
      end if
!
      ta(1) = 0.0_real64
      do i = 2, mtp
         ta(i) = (radial_large_components(i,ia)*radial_large_components(i,ic) + radial_small_components(i,ia)* &
            & radial_small_components(i,ic))*radial_derivative_over_radius(i)*tb(i)
      end do
!
      call quad (result)
      slater = result
!
!   Debug printout
!
      if (ldbpr(10)) write (99, 300) k, orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
         & orbital_symmetry_labels(ib), orbital_principal_n(ic)&
         , orbital_symmetry_labels(ic), orbital_principal_n(id), orbital_symmetry_labels(id), slater
!
      return
!
  300 format(/,"  (",1i1,")",/," R   (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2&
         ,1a2,") ","= ",1p,d19.12)
      return
!
      end function slater
end module slater_i
