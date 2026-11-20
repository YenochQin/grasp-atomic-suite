!***********************************************************************
!                                                                      *
      subroutine setout_gj(name, nci, gjc_diag, dgjc_diag)
!                                                                      *
!   Writes diagonal Landé-factor data to <name>.gj or <name>.cgj.     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state,           only: vacuum_speed_of_light_au
      use angular_symmetry_labels,         only: labj=>jlbr, labp=>jlbp
      use isotope_properties,        only: sqn, dmomnm, qmomb,                  &
                                 hfsi=>sqn, hfsd=>dmomnm, hfsq=>qmomb
      use selected_level_indices,          only: nvec, ivec
      use level_symmetry_metadata,          only: iatjpo, iaspar
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_open, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=24), intent(in) :: name
      integer, intent(in) :: nci
      real(real64), intent(in) :: gjc_diag(nvec), dgjc_diag(nvec)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ierr, jj, k
      real(real64) :: fj, gja1, gj, dgj
      character :: filnam*256, form*11, status*3
!-----------------------------------------------
!
      k = index(name," ")
      if (nci == 0) then
         filnam = name(1:k-1)//".cgj"
      else
         filnam = name(1:k-1)//".gj"
      end if
!
      form = "FORMATTED"
      status = "NEW"
      call openfl(29, filnam, form, status, ierr)
      if (ierr /= 0) then
         write (6, *) "Error when opening", filnam
         stop
      end if
!
      write (29, '("Nuclear spin                        ",1p,d22.15," au")') hfsi
      write (29, '("Nuclear magnetic dipole moment      ",1p,d22.15," n.m.")') hfsd
      write (29, '("Nuclear electric quadrupole moment  ",1p,d22.15," barns")') hfsq
      write (29, '(/,/," Interaction constants:"/,/," Level1  J Parity ",10x,"g_J",14x,"delta g_J",11x,"total g_J"/&
         &)')
!
      do i = 1, nvec
         jj = iatjpo(i)
         if (jj <= 1) cycle
         fj = 0.5_real64*dble(jj - 1)
         gja1 = sqrt(1.0_real64/(fj*(fj + 1.0_real64)))
         gj = vacuum_speed_of_light_au*gja1*gjc_diag(i)
         dgj = 0.001160_real64*gja1*dgjc_diag(i)
!
         write (29, '(1x,1i3,5x,2a4,1p,3d20.10)') ivec(i), labj(jj), labp((iaspar(i)+3)/2), gj, dgj, &
            gj + dgj
      end do
!
      close(29)
      return
!
      end subroutine setout_gj
