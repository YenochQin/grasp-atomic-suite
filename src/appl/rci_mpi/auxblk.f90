!************************************************************************
     subroutine auxblk(j2max, atwinv)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals
      use decide_c
      use def_c
      use grid_c
      use ncdist_c
      use tatb_c
      use coeils_c
      use vpilst_c,indoei=>indvpi, valoei=>valvpi
      use bilst_c
      use keilst_c
      use vinlst_c
      use orb_c, only: csf_count, orbital_count, csf_subshell_occupations
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ncharg_i
      use vacpol_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j2max
      real(real64), intent(out) :: atwinv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k
      logical :: ldbpg
!-----------------------------------------------
!
      frstco = .true.
      ncoei = 0

      if (ltrans) then
!        ...Check the maximum numbers of orbtitals allowed in brint.f
         i = max_orbitals
         select case (j2max)
         case (11)
            i = 114
         case (12)
            i = 112
         case (13)
            i = 110
         case (14)
            i = 108
         case (15)
            i = 106
         case (16)
            i = 105
         case (17)
            i = 103
         case (18)
            i = 101
         case (19)
            i = 100
         case (21:)
            i = 90
         end select

         if (i < orbital_count) then
            write (stderr_unit, *) 'In setham. The number of orbitals is too'
            write (stderr_unit, *) 'large for the brint routine'
            stop
         endif

         first = .true.
         ntpi = 0
      endif
!
!   Initialisations for the vacuum polarisation corrections
!
      if (lvp) then
         call ncharg
         call vacpol
         zdist(2:radial_point_count) = tb(2:radial_point_count)*radial_grid_derivative(2:radial_point_count)
         frstvp = .true.
         nvpi = 0
      endif
!
!   Initialisations for nuclear translational energy corrections
!
      if (emn > 0._real64) then
         atwinv = 1._real64/emn
         if (lnms) then
            frstki = .true.
            nkei = 0
         endif
         if (lsms) then
            frstvi = .true.
            nvinti = 0
         endif
      else
         ! atwinv will not be used
         lnms = .false.
         lsms = .false.
      endif

      return
      end subroutine auxblk
