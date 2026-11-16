!***********************************************************************
      subroutine orbout(rwffile2)
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grid_c
      use orb_c
      use wave_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: rwffile2*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, mfj, i
!-----------------------------------------------
!
      open(23, file=rwffile2, status="UNKNOWN", form="UNFORMATTED", position=&
         "asis")
      write (23) "G92RWF"
      do j = 1, orbital_count
         mfj = orbital_last_grid_index(j)
         write (23) orbital_principal_n(j), orbital_kappa(j), orbital_energy_parameters(j), mfj
         write (23) orbital_origin_parameters(j), (radial_large_components(i,j),i=1,mfj), (radial_small_components(i,j),i=1,mfj)
         write (23) (radial_grid(i),i=1,mfj)               ! This is a waste of resources
      end do
      close(23)

      return
      end subroutine orbout
