!***********************************************************************
!                                                                      *
      subroutine intrpqf(source_large_component, source_small_component, source_grid_size, source_radial_grid, orbital_index, &
         & original_norm_squared)
!                                                                      *
!   This  subprogram  interpolates  the  arrays  PA(1:MA), QA(1:MA),   *
!   tabulated on grid RA(1:MA) into the COMMON arrays PF(1:MF(J),J),   *
!   QF(1:MF(J),J). (Aitken's  algorithm is used. See F B Hildebrand,   *
!   Introduction  to  Numerical  Analysis, 2nd ed., McGraw-Hill, New   *
!   York, NY, 1974.) The orbital is renormalized.                      *
!                                                                      *
!   SUBROUTINEs called: RINT.                                          *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use orbital_csf_basis, only: max_orbitals
      use atomic_calculation_state,           only:accy
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rintff_i, only: rintff
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: source_grid_size
      integer, intent(in) :: orbital_index
      real(real64), intent(out) :: original_norm_squared
      real(real64), intent(in) :: source_large_component(*)
      real(real64), intent(in) :: source_small_component(*)
      real(real64), intent(in) :: source_radial_grid(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: mxord = 13
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, last_radial_index, nrstlo, kount, irow, k, nrsthi, llo, lhi, locnxt, &
         ilirok, ildiag, ilothr, mfjp1
      real(real64) :: x(mxord), dx(mxord)
      real(real64) :: polyp((mxord*(mxord + 1))/2), polyq((mxord*(mxord + 1))/2)
      real(real64) :: rama, rn, xbar, pestl, qestl, diff, difft, dxkmn1, dxirow&
         , factor, pestt, qestt, dpbp, dqbq, dnfac
      logical :: set
      logical :: used(max_radial_points)
!-----------------------------------------------
!
!   MXORD is the maximum order of the interpolation
!
!
!   This function dispenses with the need for a two-dimensional
!   array
!
!
!   Initialization
!
      rama = source_radial_grid(source_grid_size)
      rn = radial_grid(radial_point_count)
!
!   This is always true in GRASP
!
      pfff(1,orbital_index) = 0.0_real64
      qfff(1,orbital_index) = 0.0_real64
!
!   Checks
!
      if (rama > rn) then
         write (*, '(/,"INTRPQ: Grid of insufficient extent:"/," Present grid has R(N) = ",1p,1d19.12,&
            &" Bohr radii"/,"          Require R(N) = ",1d19.12," Bohr radii")') rn, rama
         stop
      end if
!
!   Determine end of grid
!
      i = radial_point_count
    1 continue
      i = i - 1
      if (radial_grid(i) <= rama) then
         last_radial_index = i
      else
         goto 1
      end if
      mfff(orbital_index) = last_radial_index
!
!   Overall initialization for interpolation
!
      nrstlo = 0
      kount = 0
!
!   Perform interpolation
!
      do i = 2, last_radial_index
!
!   Initialization for interpolation
!
         xbar = radial_grid(i)
         irow = 0
         pestl = 0.0_real64
         qestl = 0.0_real64
!
!   Determine the nearest two grid points bounding the present
!   grid point
!
    2    continue
         k = nrstlo + 1
         if (source_radial_grid(k) < xbar) then
            nrstlo = k
            goto 2
         else
            nrsthi = k
         end if
!
!   Clear relevant piece of use-indicator array
!
         llo = max(nrstlo - mxord,1)
         lhi = min(nrsthi + mxord,source_grid_size)
         used(llo:lhi) = .false.
!
!   Determine next-nearest grid point
!
    4    continue
         irow = irow + 1
         llo = max(nrstlo - irow + 1,1)
         lhi = min(nrsthi + irow - 1,source_grid_size)
         set = .false.
         do k = llo, lhi
            if (used(k)) cycle
            if (.not.set) then
               diff = source_radial_grid(k) - xbar
               locnxt = k
               set = .true.
            else
               difft = source_radial_grid(k) - xbar
               if (abs(difft) < abs(diff)) then
                  diff = difft
                  locnxt = k
               end if
            end if
         end do
         used(locnxt) = .true.
         x(irow) = source_radial_grid(locnxt)
         dx(irow) = diff
!
!   Fill table for this row
!
         do k = 1, irow
            ilirok = iloc(irow,k)
            if (k == 1) then
               polyp(ilirok) = source_large_component(locnxt)
               polyq(ilirok) = source_small_component(locnxt)
            else
               ildiag = iloc(k - 1,k - 1)
               ilothr = iloc(irow,k - 1)
               dxkmn1 = dx(k-1)
               dxirow = dx(irow)
               factor = 1.0_real64/(x(irow)-x(k-1))
               polyp(ilirok) = (polyp(ildiag)*dxirow-polyp(ilothr)*dxkmn1)*&
                  factor
               polyq(ilirok) = (polyq(ildiag)*dxirow-polyq(ilothr)*dxkmn1)*&
                  factor
            end if
         end do
!
!   Check for convergence
!
         ildiag = iloc(irow,irow)
         pestt = polyp(ildiag)
         qestt = polyq(ildiag)
         if (pestt==0.0_real64 .or. qestt==0.0_real64) then
            if (irow < mxord) then
               goto 4
            else
               pfff(i,orbital_index) = pestt
               qfff(i,orbital_index) = qestt
            end if
         else
            dpbp = abs((pestt - pestl)/pestt)
            dqbq = abs((qestt - qestl)/qestt)
            if (dqbq<accy .and. dpbp<accy) then
               pfff(i,orbital_index) = pestt
               qfff(i,orbital_index) = qestt
            else
               pestl = pestt
               qestl = qestt
               if (irow < mxord) then
                  goto 4
               else
                  pfff(i,orbital_index) = pestt
                  qfff(i,orbital_index) = qestt
                  kount = kount + 1
               end if
            end if
         end if
!
      end do
!
!   Ensure that all points of the array are defined by setting the
!   tail to zero
!
      mfjp1 = last_radial_index + 1
      pfff(mfjp1:radial_point_count,orbital_index) = 0.0_real64
      qfff(mfjp1:radial_point_count,orbital_index) = 0.0_real64
!
      if (ldbpr(3) .and. kount>0) write (99, '(/,"INTRPQ: Interpolation procedure not converged to",1p,1d19.12,&
         &" for ",1i3," of ",1i3," tabulation points")') accy, kount, last_radial_index
!
!   Normalization
!
      original_norm_squared = rintff(orbital_index,orbital_index,0)
      dnfac = 1.0_real64/sqrt(original_norm_squared)
      pfff(:last_radial_index,orbital_index) = pfff(:last_radial_index,orbital_index)*dnfac
      qfff(:last_radial_index,orbital_index) = qfff(:last_radial_index,orbital_index)*dnfac
!
      return
!
      return
      contains


      integer function iloc (ind1, ind2)
      integer, intent(in) :: ind1
      integer, intent(in) :: ind2
      iloc = (ind1*(ind1 - 1))/2 + ind2
      return
      end function iloc
!
      end subroutine intrpqf
