!***********************************************************************
!                                                                      *
      program gj
!                                                                      *
!   Entry routine for the serial Landé-factor program.                 *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use default_c
      use iounit_c
      use memory_man
      use prnt_c,          only: nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use compute_gj_only_i
      use gethfd_i
      use getmixblock_i
      use getyn_i
      use setcon_i
      use setcsla_i
      use setdbg_i
      use setmc_i
      use setout_gj_i
      use factt_i
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, nci, ncore_not_used
      logical :: yes
      real(real64), pointer :: gjc_diag(:), dgjc_diag(:)
      character :: name*24
!-----------------------------------------------
!
      write (stderr_unit, *)
      write (stderr_unit, *) 'GJ'
      write (stderr_unit, *) 'This is the serial Lande-factor program'
      write (stderr_unit, *) 'Input files:  isodata, name.c, name.(c)m, name.w'
      write (stderr_unit, *) 'Output files: name.gj or name.cgj'
!
      write (stderr_unit, *)
      write (stderr_unit, *) 'Default settings?'
      yes = getyn()
      write (stderr_unit, *)
      if (yes) then
         ndef = 0
      else
         ndef = 1
      endif
!
   10 continue
      write (stderr_unit, *) 'Name of state'
      read (*, '(A)') name
      k = index(name,' ')
      if (k == 1) then
         write (stderr_unit, *) 'Names may not start with a blank'
         go to 10
      endif
!
      write (stderr_unit, *)
      write (stderr_unit, *) 'Mixing coefficients from a CI calc.?'
      yes = getyn()
      if (yes) then
         nci = 0
      else
         nci = 1
      endif
!
      call setdbg
      call setmc
      call setcon
      call setcsla(name, ncore_not_used)
      call gethfd(name)
      call getmixblock(name, nci)
      call factt
!
      call alloc(gjc_diag, nvec, 'GJC_DIAG', 'GJ')
      call alloc(dgjc_diag, nvec, 'DGJC_DIAG', 'GJ')
!
      call compute_gj_only(gjc_diag, dgjc_diag)
      call setout_gj(name, nci, gjc_diag, dgjc_diag)
!
      call dalloc(gjc_diag, 'GJC_DIAG', 'GJ')
      call dalloc(dgjc_diag, 'DGJC_DIAG', 'GJ')
!
      write (stderr_unit, *)
      write (stderr_unit, *) 'GJ: Execution complete.'
      stop
      end program gj
