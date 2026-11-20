!************************************************************************
      subroutine cpmix(name, inpci)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use array_memory_management, only: alloc, dalloc
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
      integer :: copied_mixing_unit
      integer :: original_mixing_unit
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: inpci
      character, intent(inout) :: name(2)*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: io_status, nelec, total_csf_count, nw, nvectot, vector_size, nblock, iblk, ib, &
         ncf, nvec, iatjp, iaspa, i, j
      real(real64), pointer :: eval(:), evec(:)
      integer, pointer :: ivec(:)
      real(real64) :: hamiltonian_energy_shift
      character :: mixing_file_signature*6
!-----------------------------------------------



      name(2) = trim(name(1))//"_CP"
      if (inpci == 0) then
         open(newunit=copied_mixing_unit, file=trim(name(2))//".cbm", form="UNFORMATTED", status=&
            "UNKNOWN", position="asis")
      else
         open(newunit=copied_mixing_unit, file=trim(name(2))//".bm", form="UNFORMATTED", status=&
            "UNKNOWN", position="asis")
      end if
      if (inpci == 0) then
         open(newunit=original_mixing_unit, file=trim(name(1))//".cbm", form="UNFORMATTED", status=&
            "OLD", position="asis")
      else
         open(newunit=original_mixing_unit, file=trim(name(1))//".bm", form="UNFORMATTED", status=&
            "OLD", position="asis")
      end if
      read (original_mixing_unit, iostat=io_status) mixing_file_signature
      if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
         write (*, *) "Not a GRASP mixing file"
         stop
      end if
      write (copied_mixing_unit) mixing_file_signature
      read (original_mixing_unit) nelec, total_csf_count, nw, nvectot, vector_size, nblock
      write (copied_mixing_unit) nelec, total_csf_count, nw, nvectot, vector_size, nblock

      do iblk = 1, nblock
         read (original_mixing_unit) ib, ncf, nvec, iatjp, iaspa
         write (copied_mixing_unit) ib, ncf, nvec, iatjp, iaspa
         call alloc (eval, nvec, "EVAL", "CPMIX")
         call alloc (evec, ncf*nvec, "EVEC", "CPMIX")
         call alloc (ivec, nvec, "IVEC", "CPMIX")
         read (original_mixing_unit) (ivec(i),i=1,nvec)
         read (original_mixing_unit) hamiltonian_energy_shift, (eval(i),i=1,nvec)
         read (original_mixing_unit) ((evec(i + (j - 1)*ncf),i=1,ncf),j=1,nvec)

         write (copied_mixing_unit) (ivec(i),i=1,nvec)
         write (copied_mixing_unit) hamiltonian_energy_shift, (eval(i),i=1,nvec)
         write (copied_mixing_unit) ((evec(i + (j - 1)*ncf),i=1,ncf),j=1,nvec)
         call dalloc (eval, "EVAL", "CPMIX")
         call dalloc (evec, "EVEC", "CPMIX")
         call dalloc (ivec, "IVEC", "CPMIX")
      end do
      name(2) = name(1)
      close(original_mixing_unit)
      close(copied_mixing_unit)
      return
      end subroutine cpmix
