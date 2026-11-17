!***********************************************************************
!                                                                      *
      subroutine radfile(name)
!                                                                      *
!  This subroutine outputs the transformed radial orbitals             *
!                                                                      *
!  Written by Per Jonsson                                              *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use biorb_c, only: double, max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use grid_c, only: double, radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use orb_c, only: double, max_orbitals
      use wave_c, only: double, max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character, intent(in) :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, k, i
!-----------------------------------------------
!
      j = index(name(1)," ")
      open(unit=30, file=name(1)(1:j-1)//".bw", form="UNFORMATTED", status=&
         "UNKNOWN", position="asis")

      write (30) "G92RWF"
      write (*, *) "NWII", nwii
      do k = 1, nwii
         write (30) npii(k), nakii(k), eii(k), mfii(k)
         write (30) pzii(k), (pfii(i,k),i=1,mfii(k)), (qfii(i,k),i=1,mfii(k))
         write (30) (radial_grid(i),i=1,mfii(k))
      end do

      close(30)

      j = index(name(2)," ")
      open(unit=30, file=name(2)(1:j-1)//".bw", form="UNFORMATTED", status=&
         "UNKNOWN", position="asis")

      write (30) "G92RWF"
      write (*, *) "NWFF", nwff
      do k = 1, nwff
         write (30) npff(k), nakff(k), eff(k), mfff(k)
         write (30) pzff(k),(pfff(i,k),i=1,mfff(k)),(qfff(i,k),i=1,mfff(k))
         write (30) (radial_grid(i),i=1,mfff(k))
      end do

      close(30)

      return
      end subroutine radfile
