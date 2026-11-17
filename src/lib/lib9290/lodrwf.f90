module lodrwf_i
    implicit none (type, external)
    private
    public :: lodrwf
contains
!***********************************************************************
!                                                                      *
      subroutine lodrwf(iierr)
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
!   Block version by Xinghong He          Last revision: 27 May 1997   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      use memory_man
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use def_c,           only: speed_of_light_au, nuclear_charge
      use grid_c
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use orb_c
      use wave_c,          only: orbital_origin_parameters, radial_large_components,radial_small_components
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use intrpq_i, only: intrpq
      use orthsc_i, only: orthsc
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: iierr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, k, nwin, ios, npy, naky, my
      integer :: ierr
      real(real64) :: con, fkk, ey, pzy, dnorm
      real(real64), pointer :: pa(:), qa(:), ra(:)
!-----------------------------------------------
!
!
!   Write entry message
!
      write (6, *) "Loading Radial WaveFunction File ..."
!
!   Allocate storage to orbital arrays
!
      call alloc(radial_large_components, max_radial_points, orbital_count, "PF", "LODRWF")
      call alloc(radial_small_components, max_radial_points, orbital_count, "QF", "LODRWF")

!
!   Setup: (1) Orbital arrays to zero
!          (2) Array E to -1 (no orbitals estimated)
!          (3) Parameters GAMMA for each orbital
!
      con = nuclear_charge/speed_of_light_au
      con = con*con
!
      do j = 1, orbital_count
!
         radial_large_components(:radial_point_count,j) = 0.0_real64
         radial_small_components(:radial_point_count,j) = 0.0_real64
!
         orbital_energy_parameters(j) = -1.0_real64
!
         k = abs(orbital_kappa(j))
         if (nuclear_model_parameter_count > 0) then
            orbital_origin_exponents(j) = dble(k)
         else if (nuclear_model_parameter_count == 0) then
            fkk = dble(k*k)
            if (fkk >= con) then
               orbital_origin_exponents(j) = sqrt(fkk - con)
            else
               !WRITE (istde,*) 'LODRWF: Imaginary gamma parameter'
               !WRITE (istde,*) ' for ',NP(J),NH(J),' orbital; the'
               !WRITE (istde,*) ' point model for the nucleus'
               !WRITE (istde,*) ' is inappropriate for Z > ',C,'.'
               stop "lodrwf: Inappropriate gamma"
            end if
         end if
!
      end do
!
!   Read orbital information from Read Orbitals File; write summary
!   to  .dbg  file if option set
!
      if (ldbpr(3)) write (99, 300)
      nwin = 0
    3 continue
      read (23, iostat=ios) npy, naky, ey, my

      if (ios == 0) then
         call alloc (pa,my, "PA", "LODRWF")
         call alloc (qa,my, "QA", "LODRWF")
         call alloc (ra,my, "RA", "LODRWF")

         read (23) pzy, (pa(i),i=1,my), (qa(i),i=1,my)
         read (23) (ra(i),i=1,my)

         do j = 1, orbital_count
            if (.not.(orbital_energy_parameters(j)<0.0_real64 .and. npy==orbital_principal_n(j) .and. naky==orbital_kappa(j))) cycle
            orbital_origin_parameters(j) = pzy
            orbital_energy_parameters(j) = ey
            call intrpq(pa, qa, my, ra, j, dnorm)
            if (ldbpr(3)) write (99, 301) orbital_principal_n(j), orbital_symmetry_labels(j), orbital_energy_parameters(j), dnorm
            nwin = nwin + 1
         end do
         call dalloc (pa, "PA", "LODRWF" )
         call dalloc (qa, "QA", "LODRWF")
         call dalloc (ra, "RA", "LODRWF")
         goto 3
      end if
      if (ldbpr(3)) write (99, *) " orbitals renormalised;"
!
!   Stop with an error message if all orbitals are not known
!
      if (nwin < orbital_count) then
         iierr = 1
         goto 5
      end if
!
!   Schmidt orthogonalise the orbitals
!
      call orthsc
      if (ldbpr(3)) write (99, *) " orbitals orthogonalised and renormalised;"
!
      iierr = 0
    5 continue
      return
!
  300 format(/,"From SUBROUTINE LODRWF:"/," Orbital",8x,"Eigenvalue",19x,"Norm")
  301 format(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)
      return
!
      end subroutine lodrwf
end module lodrwf_i
