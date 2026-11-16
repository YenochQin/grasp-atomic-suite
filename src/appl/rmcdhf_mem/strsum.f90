!***********************************************************************
      subroutine strsum

!   Generates the first part of  rscf92.sum  (on stream 24).
!
!   Call(s) to: [LIB92] CALEN, CONVRT.
!
!   Written by Farid A. Parpia            Last revision: 26 Sep 1993
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
      use foparm_c
      use grid_c
      use mcpb_c
      use npar_c
      use npot_c
      use orb_c
      use wave_c
      use wfac_c
!      USE iccu_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: lenth, iend, i, ibeg
      character :: record*256, cdata*26, ctime*8, cdate*8, clevel*4
!-----------------------------------------------
!
!     POINTER (PWEIGH,WEIGHT(1))
!     POINTER (PCCMIN,ICCMIN(1))
!     POINTER (PNTRPF,PF(NNNP,1))
!     POINTER (PNTRQF,QF(NNNP,1))
!
!   Both the nuclear charge and the number of electrons are
!   known at this point; load IONCTY with the ionicity
!
      ioncty = nint(nuclear_charge) - electron_count
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!      CALL CALEN (CTIME, CDATE)
!      WRITE (24,*) 'RSCF92 run at ',CTIME,' on ',CDATE,'.'
!
!   Write out the basic dimensions of the electron cloud
!
      write (24, *)
      call convrt (electron_count, record, lenth)
      write (24, *) 'There are '//record(1:lenth)//' electrons in the cloud'
      call convrt (csf_count, record, lenth)
      write (24, *) ' in '//record(1:lenth)//' relativistic CSFs'
      call convrt (orbital_count, record, lenth)
      write (24, *) ' based on '//record(1:lenth)//' relativistic subshells.'
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      if (lfordr) then
         write (24, *)
         call convrt (iccut, record, lenth)
         write (24, *) ' CSFs 1--'//record(1:lenth)//' constitute'//&
            ' the zero-order space;'
      endif
!
!   Write out the nuclear parameters
!
      write (24, *)
      write (24, 300) nuclear_charge

      if (emn == 0._real64) then
         write (24, *) ' the nucleus is stationary;'
      else
         write (24, 301) emn
      endif
      if (nuclear_model_parameter_count == 2) then
         write (24, *) ' Fermi nucleus:'
         write (24, 302) nuclear_model_parameters(1), nuclear_model_parameters(2)
         call convrt (first_grid_index_beyond_fermi_c, record, lenth)
         write (24, *) ' there are '//record(1:lenth)//&
            ' tabulation points in the nucleus.'
      else
         write (24, *) ' point nucleus.'
      endif
!
!   Write out the physical effects specifications
!
      write (24, *)
      write (24, 303) speed_of_light_au
!
!   Write out the parameters of the radial grid
!
      write (24, *)
      if (asymptotic_grid_step == 0._real64) then
         write (24, 305) radial_grid_scale, transformed_grid_step, radial_point_count
      else
         write (24, 306) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
      endif
      write (24, 307) radial_grid(1), radial_grid(2), radial_grid(radial_point_count)
      write (24, *)
!
!  (E)AL calculation, returns here
!
      if (ncmin == 0) then
         write (24, *) '(E)AL calculation.'
         return
      endif
!
!  Info exclusively for EOL calculations
!
      if (ncmin == 1) then
         write (24, *) 'OL calculation.'
         call convrt (iccmin(1), record, lenth)
         write (24, *) 'Level '//record(1:lenth)//' will be optimised.'
      else
         write (24, *) 'EOL calculation.'
         call convrt (ncmin, record, lenth)
         write (24, *) record(1:lenth)//' levels will be optimised;'
         record(1:20) = ' their indices are: '
         iend = 20
         do i = 1, ncmin
            ibeg = iend + 1
            call convrt (iccmin(i), clevel, lenth)
            if (i /= ncmin) then
               iend = ibeg + lenth + 1
               record(ibeg:iend) = clevel(1:lenth)//', '
            else
               iend = ibeg + lenth
               record(ibeg:iend) = clevel(1:lenth)//'.'
            endif
            if (iend < 120) cycle
            write (24, *) record(1:iend)
            record(1:2) = '  '
            iend = 2
         end do
         if (iend /= 2) write (24, *) record(1:iend)
         if (weight(1) == (-1._real64)) then
            write (24, *) 'Each is assigned its statistical weight;'
         else if (weight(1) == (-2._real64)) then
            write (24, *) 'All levels are weighted equally;'
         else
            write (24, *) ' weighted as follows:'
            write (24, *) (weight(i),i=1,ncmin)
         endif
      endif

  300 format('The atomic number is ',1f14.10,';')
  301 format(' the mass of the nucleus is ',1p,d19.12,' electron masses;')
  302 format('  c =',1p,1d19.12,' Bohr radii,'/,'  a =',1d19.12,' Bohr radii;')
  303 format('Speed of light = ',3p,d19.12,' atomic units.')
  305 format('Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',' I = 1, ..., N;'/,/,&
         ' RNT  = ',1p,d19.12,' Bohr radii;'/,' H    = ',d19.12,' Bohr radii;'/&
         ,' N    = ',1i4,';')
  306 format('Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',&
         ' I = 1, ..., N;'/,/,' RNT  = ',1p,d19.12,' Bohr radii;'/,' H    = ',d&
         19.12,' Bohr radii;'/,' HP   = ',d19.12,' Bohr radii;'/,' N    = ',1i4&
         ,';')
  307 format(' R(1) = ',1p,1d19.12,' Bohr radii;'/,' R(2) = ',1d19.12,&
         ' Bohr radii;'/,' R(N) = ',1d19.12,' Bohr radii.')

      return
      end subroutine strsum
