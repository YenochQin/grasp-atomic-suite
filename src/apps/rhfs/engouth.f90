!***********************************************************************
!                                                                      *
      subroutine engouth(hamiltonian_energy_shift, e, jtot, ipar, level_indices, level_count, mode)
      use calculation_file_units, only: hyperfine_output_unit, summary_output_unit
!                                                                      *
!   This  subroutine prints  energy levels, splittings, and energies   *
!   relative to the lowest in  Hartrees, Kaysers, and  eV, using the   *
!   reduced mass corrected value for the Rydberg. If MODE is 0, only   *
!   the eigenenergies are printed. If  MODE  is 1, the eigenenergies   *
!   and separations are printed. If MODE is 2, the eigenenergies and   *
!   energies relative to level 1 are printed. If MODE is 3, the eig-   *
!   enenergies,  separations,  and energies relative to level  1 are   *
!   printed.                                                           *
!                                          Last updated: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: &
          hartree_to_wavenumber_cm, hartree_to_ev, speed_of_light_cm_per_s, einstein_a_au_to_si, einstein_b_au_to_si
      use angular_symmetry_labels, only: right_justified_j_labels => jlbr, parity_labels => jlbp
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: level_count
      integer, intent(in) :: mode
      real(real64), intent(in) :: hamiltonian_energy_shift
      integer, intent(in) :: jtot(level_count)
      integer, intent(in) :: ipar(level_count)
      integer, intent(in) :: level_indices(level_count)
      real(real64), intent(in) :: e(level_count)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: j, i, ip
      real(real64) :: eau, ecm, eev
!-----------------------------------------------
!
!   Always print the eigenenergies
!
      write (summary_output_unit, '(/,"Eigenenergies:")')
      write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
      write (hyperfine_output_unit, '(/,"Eigenenergies:")')
      write (hyperfine_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
      do j = 1, level_count
         i = level_indices(j)
         eau = e(j) + hamiltonian_energy_shift
         ecm = eau*hartree_to_wavenumber_cm
         eev = eau*hartree_to_ev
         ip = (ipar(j)+3)/2
         write (summary_output_unit, '(1i3,2x,2a4,1p,3d22.14)') i, right_justified_j_labels(jtot(j)), parity_labels(ip), eau, &
            & ecm, eev
         write (hyperfine_output_unit, '(1i3,2x,2a4,1p,3d22.14)') i, right_justified_j_labels(jtot(j)), parity_labels(ip), eau, &
            & ecm, eev
      end do
!
      if (level_count > 1) then
!
!   Energy separations
!
         if (mode==1 .or. mode==3) then
            write (summary_output_unit, '(/,"Energy of each level relative to immediately lower"," level:")')
            write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
            write (hyperfine_output_unit, '(/,"Energy of each level relative to immediately lower"," level:")')
            write (hyperfine_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
            do j = 2, level_count
               i = level_indices(j)
               eau = e(j) - e(j-1)
               ecm = eau*hartree_to_wavenumber_cm
               eev = eau*hartree_to_ev
               ip = (ipar(j)+3)/2
               write (summary_output_unit, '(1i3,2x,2a4,1p,3d22.14)') i, right_justified_j_labels(jtot(j)), parity_labels(ip), &
                  & eau, ecm, eev
               write (hyperfine_output_unit, '(1i3,2x,2a4,1p,3d22.14)') i, right_justified_j_labels(jtot(j)), &
                  & parity_labels(ip), eau, ecm, eev
            end do
         end if
!
!   Energies relative to level 1
!
         if (mode==2 .or. mode==3) then
            write (summary_output_unit, '(/,"Energy of each level relative to lowest level:")')
            write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
            write (hyperfine_output_unit, '(/,"Energy of each level relative to lowest level:")')
            write (hyperfine_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
            do j = 2, level_count
               i = level_indices(j)
               eau = e(j) - e(1)
               ecm = eau*hartree_to_wavenumber_cm
               eev = eau*hartree_to_ev
               ip = (ipar(j)+3)/2
               write (summary_output_unit, '(1i3,2x,2a4,1p,3d22.14)') i, right_justified_j_labels(jtot(j)), parity_labels(ip), &
                  & eau, ecm, eev
               write (hyperfine_output_unit, '(1i3,2x,2a4,1p,3d22.14)') i, right_justified_j_labels(jtot(j)), &
                  & parity_labels(ip), eau, ecm, eev
            end do
         end if
!
      end if
!
      return
!
      end subroutine engouth
