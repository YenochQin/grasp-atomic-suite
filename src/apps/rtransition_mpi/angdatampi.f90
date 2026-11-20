!***********************************************************************
!                                                                      *
      subroutine angdata(name, avail, jkp, transition_angular_unit)
      use, intrinsic :: iso_fortran_env, only: output_unit
!                                                                      *
!   Checks if the angular file name(1).name(2).T is available          *
!   and appropriate                                                    *
!                                                                      *
!   Written by Per Jonsson                      6 March 1997           *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations
      use transition_matrix_workspace, only: nkp, kp
      use mpi_state, only: mpi_process_count, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: jkp
      integer, intent(out) :: transition_angular_unit
      logical, intent(out) :: avail
      character, intent(inout) :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j1,j2,iblki,iblkf,nwd,nkpd,nprocsd,myidd
      logical :: found
      character :: s(-9:9)*2
!-----------------------------------------------

      s((-9)) = "-9"
      s((-8)) = "-8"
      s((-7)) = "-7"
      s((-6)) = "-6"
      s((-5)) = "-5"
      s((-4)) = "-4"
      s((-3)) = "-3"
      s((-2)) = "-2"
      s((-1)) = "-1"
      s(0) = "+0"
      s(1) = "+1"
      s(2) = "+2"
      s(3) = "+3"
      s(4) = "+4"
      s(5) = "+5"
      s(6) = "+6"
      s(7) = "+7"
      s(8) = "+8"
      s(9) = "+9"

      j1 = index(name(1)," ")
      j2 = index(name(2)," ")
      inquire(file=name(1)(1:j1-1)//"."//name(2)(1:j2-1)//"."//s(kp(jkp))//"T"&
         , exist=found)
      if (.not.found) then
         write (output_unit, *)
         write (output_unit, *) " Angular file not available"
         avail = .false.
         return
      else
!
!  Open the file and check if it is appropriate for the present case
!
         open(newunit=transition_angular_unit, file=name(1)(1:j1-1)//"."//name(2)(1:j2-1)//"."//s(&
            kp(jkp))//"T", status="OLD", form="UNFORMATTED", position="asis")
         rewind (transition_angular_unit)
         read (transition_angular_unit) iblki,iblkf,nwd,nkpd,nprocsd,myidd
         if (.not.(nwd==orbital_count .and. nkpd==nkp                          &
             .and.nprocsd==mpi_process_count.and.myidd==mpi_rank)) then
            write (output_unit, *) " Angular file not appropriate"
            avail = .false.
            close(transition_angular_unit, status="DELETE")
            return
         else
            rewind (transition_angular_unit)
            write (output_unit, *) " Angular data read from file"
            avail = .true.
         end if
      end if
      return
      end subroutine angdata
