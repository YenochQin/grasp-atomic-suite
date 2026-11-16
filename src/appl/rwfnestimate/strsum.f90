

!***********************************************************************
!                                                                      *
      subroutine strsum
!                                                                      *
!   Generates the first part of  erwf.sum  (on stream 24).             *
!                                                                      *
!   Call(s) to: [LIB92] CALEN, CONVRT.                                 *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
      use grid_c
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use npot_c, only: first_grid_index_beyond_fermi_c
      use orb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use calen_i, only: calen
      use convrt_i, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: lenth
      character(len=10) :: ctime
      character(len=8) :: cdate
      character(len=256) :: record
      character(len=26)  :: cdata

!-----------------------------------------------
!
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
      call calen (ctime, cdate)
      write (24, *) "ERWF run at ", ctime, " on ", cdate, "."
!
!   Write out the basic dimensions of the electron cloud
!
      write (24, *)
      call convrt (electron_count, record, lenth)
      write (24, *) "There are "//record(1:lenth)//" electrons in the cloud"
      call convrt (csf_count, record, lenth)
      write (24, *) " in "//record(1:lenth)//" relativistic CSFs"
      call convrt (orbital_count, record, lenth)
      write (24, *) " based on "//record(1:lenth)//" relativistic subshells."
!
!   Write out the nuclear parameters
!
      write (24, *)
      write (24, 300) nuclear_charge
      if (nuclear_model_parameter_count == 2) then
         write (24, *) " Fermi nucleus:"
         write (24, 302) nuclear_model_parameters(1), nuclear_model_parameters(2)
         call convrt (first_grid_index_beyond_fermi_c, record, lenth)
         write (24, *) " there are "//record(1:lenth)//&
            " tabulation points in the nucleus."
      else
         write (24, *) " point nucleus."
      end if
!
!   Write out the physical effects specifications
!
      write (24, *)
      write (24, 303) speed_of_light_au
!
!   Write out the parameters of the radial grid
!
      write (24, *)
      if (asymptotic_grid_step == 0.0_real64) then
         write (24, 305) radial_grid_scale, transformed_grid_step, radial_point_count
      else
         write (24, 306) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
      end if
      write (24, 307) radial_grid(1), radial_grid(2), radial_grid(radial_point_count)
!
      write (24, *)
!
      return
!
  300 format("The atomic number is ",1f14.10,";")
  302 format("  c =",1p,1d19.12," Bohr radii,"/,"  a =",1d19.12," Bohr radii;")
  303 format("Speed of light = ",1p,d19.12," atomic units.")
  305 format("Radial grid: R(I) = RNT*(exp((I-1)*H)-1),"," I = 1, ..., N;"/,/,&
         " RNT  = ",1p,d19.12," Bohr radii;"/," H    = ",d19.12," Bohr radii;"/&
         ," N    = ",1i4,";")
  306 format("Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,",&
         " I = 1, ..., N;"/,/," RNT  = ",1p,d19.12," Bohr radii;"/," H    = ",d&
         19.12," Bohr radii;"/," HP   = ",d19.12," Bohr radii;"/," N    = ",1i4&
         ,";")
  307 format(" R(1) = ",1p,1d19.12," Bohr radii;"/," R(2) = ",1d19.12,&
         " Bohr radii;"/," R(N) = ",1d19.12," Bohr radii.")
      return
!
      end subroutine strsum
