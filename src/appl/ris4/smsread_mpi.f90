!***********************************************************************
!                                                                      *
      subroutine smsread_mpi(vint,vint2)
!                                                                      *
!   MPI version of SMSREAD. Each rank scans the full angular file but  *
!   only accumulates records whose IC belongs to its strided subset.    *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: orbital_index_radix, max_orbitals
      use prnt_c
      use ris_c
      use orb_c
      use eigv_c
      use buffer_c
      use debug_c,          only: cutoff
      use mpi_c,            only: mpi_rank, mpi_process_count
      use alcbuf_i, only: alcbuf
      implicit none (type, external)
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      integer, parameter :: key = orbital_index_radix
      real(real64) :: evpair(nvec)
      real(real64) :: contri, contrik1, coeffsms, radial_k1, radial_sms,&
         pair_scale, term_k1, term_sms
      integer :: ic, ir, j, iia, iib, iic, iid, loc, lab, ios
      logical :: ownrec

      call alcbuf (1)
      rewind(51)
   16 read (51,iostat = ios) ic,ir
      if (ios == 0) then
         read(51) coeffsms,lab
         ownrec = mod(ic-1, mpi_process_count) == mpi_rank
         if (ownrec) then
            iid = mod (lab, key)
            lab = lab/key
            iib  = mod (lab, key)
            lab = lab/key
            iic = mod (lab, key)
            iia = lab/key
            radial_k1 = vint(iia,iic) * vint(iib,iid)
            radial_sms = (vint2(iia,iic) * vint(iib,iid)                 &
               + vint2(iib,iid) * vint(iia,iic))/2.0_real64
            if (radial_k1 == 0.0_real64 .and. radial_sms == 0.0_real64) goto 16
            term_k1 = -coeffsms * radial_k1
            term_sms = -coeffsms * radial_sms
            if (ir == ic) then
               pair_scale = 1.0_real64
            else
               pair_scale = 2.0_real64
            end if
            do j = 1,nvec
               loc = (j-1)*csf_count
               evpair(j) = pair_scale * csf_mixing_coefficients(ic+loc) * csf_mixing_coefficients(ir+loc)
            end do
            do j = 1,nvec
               contrik1 = evpair(j) * term_k1
               contri = evpair(j) * term_sms
               smsc1(j) = smsc1(j) + contrik1
               smsc2(j) = smsc2(j) + contri
            end do
         end if
         goto 16
      end if
      call alcbuf (3)
      return
      end subroutine smsread_mpi
