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
      use calculation_input_defaults, only: inpci, ndef, ndump
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use array_memory_management, only: alloc, dalloc
      use selected_level_indices,          only: nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use compute_gj_only_i, only: compute_gj_only
      use gethfd_i, only: gethfd
      use getmixblock_i, only: getmixblock
      use yes_no_input, only: getyn
      use physical_constants, only: setcon
      use csf_basis_loader, only: setcsla
      use setdbg_i, only: setdbg
      use machine_constants, only: setmc
      use setout_gj_i, only: setout_gj
      use log_factorial_table, only: factt
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
      write (stderr_unit, *) "GJ"
      write (stderr_unit, *) "This is the serial Lande-factor program"
      write (stderr_unit, *) "Input files:  isodata, name_descriptors.parquet, name_header.toml, name.(c)m, name.w"
      write (stderr_unit, *) "Output files: name.gj or name.cgj"
!
      write (stderr_unit, *)
      write (stderr_unit, *) "Default settings?"
      yes = getyn()
      write (stderr_unit, *)
      if (yes) then
         ndef = 0
      else
         ndef = 1
      end if
!
   10 continue
      write (stderr_unit, *) "Name of state"
      read (*, "(A)") name
      k = index(name," ")
      if (k == 1) then
         write (stderr_unit, *) "Names may not start with a blank"
         goto 10
      end if
!
      write (stderr_unit, *)
      write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
      yes = getyn()
      if (yes) then
         nci = 0
      else
         nci = 1
      end if
!
      call setdbg
      call setmc
      call setcon
      call setcsla(name, ncore_not_used)
      call gethfd(name)
      call getmixblock(name, nci)
      call factt
!
      call alloc(gjc_diag, nvec, "GJC_DIAG", "GJ")
      call alloc(dgjc_diag, nvec, "DGJC_DIAG", "GJ")
!
      call compute_gj_only(gjc_diag, dgjc_diag)
      call setout_gj(name, nci, gjc_diag, dgjc_diag)
!
      call dalloc(gjc_diag, "GJC_DIAG", "GJ")
      call dalloc(dgjc_diag, "DGJC_DIAG", "GJ")
!
      write (stderr_unit, *)
      write (stderr_unit, *) "GJ: Execution complete."
      stop
      end program gj
