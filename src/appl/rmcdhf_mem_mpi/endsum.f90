!***********************************************************************
!                                                                      *
      subroutine endsum
!                                                                      *
!   Generates the last part of  rscf92.sum  (on stream 24).            *
!                                                                      *
!   Call(s) to: [LIB92]: ENGOUT, RINT.                                 *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
      use eigv_c
      use orb_c
      use scf_c
!gg      USE syma_C
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rint_i, only: rint
      use engoutgg_i, only: engoutgg
      use csfwgt_i, only: csfwgt
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, mode
      real(real64) :: wa, wb, wc, wd, we
!-----------------------------------------------
!
!   Write out the orbital properties
!
      write (24, 301)
      do i = 1, orbital_count
         write (24, 302) orbital_principal_n(i),orbital_symmetry_labels(i),orbital_energy_parameters(i), &
            & orbital_origin_parameters(i),orbital_origin_exponents(i),radial_large_components(2,i),radial_small_components(2,i), &
            scnsty(i), orbital_last_grid_index(i)
      end do
!
      write (24, 303)
      do i = 1, orbital_count
         wa = rint(i,i,-1)
         wb = rint(i,i,1)
         wc = rint(i,i,2)
         wd = rint (i,i, 4)
         we = 0.0_real64
         if (orbital_symmetry_labels(i) /= "s " .and. orbital_symmetry_labels(i) /= "p-") then
            we = rint(i,i,-3)
         end if
         write (24,304) orbital_principal_n(i),orbital_symmetry_labels(i),we,wa,wb,wc,wd, ucf(i)
      end do
!
      if (ncmin /= 0) then
         mode = 0
         call engoutgg (shifted_eigenvalues,iccmin,ncmin,mode)
!GG         CALL ENGOUT (EVAL, IATJPO, IASPAR, ICCMIN, NCMIN, MODE)
         call csfwgt (.false.)
      end if
!
      close(24)
!
      return
!
  301 format(/,"Radial wavefunction summary:"/,/,67x,"Self"/,"Subshell",6x,"e",&
         13x,"p0",5x,"gamma",5x,"P(2)",7x,"Q(2)",3x,"Consistency"," MTP"/)
  302 format(1x,i2,a2,1x,1p,d17.10,1p,d11.3,0p,f6.2,1p,3(d11.3),i5)
  303 format (/18x,"-3",14x,"-1",29x,"2",14x,"4",5x,"Generalised"     &
              /"Subshell",4x,"<  r  >",8x,"<  r  >",8x,"<  r  >",8x,  &
               "<  r  >",8x,"<  r  >",6x,"occupation"/)
  304 format (1x,1i2,1a2,1x,1p,6d15.5)
      return
!
      end subroutine endsum
