!***********************************************************************
!                                                                      *
      real(real64) function brra (itype, ia, ic, ib, id, k)
!                                                                      *
!   This routine evaluates the transverse interaction integrals:       *
!                                                                      *
!      ITYPE = 1: General R (k; a c | b d )                            *
!            = 2: General S (k; a c | b d )                            *
!            = 3: R (k; a , b d )                                      *
!            = 4: F (k; a , b )                                        *
!            = 5: G (k; a , b ) type integral                          *
!            = 6: H (k; a , b ) type integral                          *
!                                                                      *
!   Call(s) to: [RCI92]: RKINT, SKINT.                                 *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      use cons_c
      use grid_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rkint_i
      use skint_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: itype
      integer, intent(in) :: ia
      integer, intent(in) :: ic
      integer, intent(in) :: ib
      integer, intent(in) :: id
      integer  :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mxrac, i, mxrbd
      real(real64) :: rac(max_radial_points), rbd(max_radial_points)
!-----------------------------------------------
!
      mxrac = min(orbital_last_grid_index(ia),orbital_last_grid_index(ic))
      do i = 1, mxrac
         rac(i) = radial_large_components(i,ia)*radial_small_components(i,ic)
      end do
!
      mxrbd = min(orbital_last_grid_index(ib),orbital_last_grid_index(id))
      do i = 1, mxrbd
         rbd(i) = radial_large_components(i,ib)*radial_small_components(i,id)
      end do
!
      go to (21,22,23,24,25,26) itype
!
!   ITYPE = 1
!
   21 continue
      if (ia==ib .and. ic==id) go to 9
      if (ia==id .and. ic==ib) go to 10
      brra = (rkint(rac,ia,ic,rbd,ib,id,k,1) + rkint(rac,ia,ic,rbd,ib,id,k,2)&
          + rkint(rbd,ib,id,rac,ia,ic,k,1) + rkint(rbd,ib,id,rac,ia,ic,k,2))*&
         half
      return
!
!   ITYPE = 2
!
   22 continue
      if (ia==ib .and. ic==id) go to 26
      if (ia==id .and. ic==ib) go to 26
      brra = (skint(rac,ia,ic,rbd,ib,id,k,1) + skint(rac,ia,ic,rbd,ib,id,k,2))*&
         half
      return
!
!   ITYPE = 3
!
   23 continue
      if (ia == ic) then
         do i = 1, mxrbd
            rbd(i) = rbd(i) + radial_large_components(i,id)*radial_small_components(i,ib)
         end do
         brra = (rkint(rac,ia,ic,rbd,ib,id,k,0) + rkint(rbd,ib,id,rac,ia,ic,k,0&
            ) + rkint(rac,ia,ic,rbd,ib,id,k,2) + rkint(rbd,ib,id,rac,ia,ic,k,2)&
            )*half
         return
      endif
      do i = 1, mxrac
         rac(i) = rac(i) + radial_large_components(i,ic)*radial_small_components(i,ia)
      end do
      brra = (rkint(rac,ia,ic,rbd,ib,id,k,1) + rkint(rbd,ib,id,rac,ia,ic,k,1)&
          + rkint(rac,ia,ic,rbd,ib,id,k,0) + rkint(rbd,ib,id,rac,ia,ic,k,0))*&
         half
      return
!
!   ITYPE = 4
!
   24 continue
      brra = rkint(rac,ia,ic,rbd,ib,id,k,0) + rkint(rbd,ib,id,rac,ia,ic,k,0)
      return
!
!   ITYPE = 5
!
   25 continue
      if (ia==id .and. ic==ib) go to 10
    9 continue
      brra = two*rkint(rac,ia,ic,rbd,ib,id,k,1)
      return
   10 continue
      brra = rkint(rac,ia,ic,rbd,ib,id,k,1) + rkint(rbd,ib,id,rac,ia,ic,k,1)
      return
!
!   ITYPE = 6
!
   26 continue
      brra = skint(rac,ia,ic,rbd,ib,id,k,1)
      return
!
      end function brra
