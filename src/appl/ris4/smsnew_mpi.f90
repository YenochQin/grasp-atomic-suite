!***********************************************************************
!                                                                      *
      subroutine smsnew_mpi(doit,vint,vint2)
!                                                                      *
!   MPI version of SMSNEW. Each rank handles a strided subset of the   *
!   outer IC loop and can optionally write a per-rank angular-data      *
!   fragment that the root process later merges into name.ITB.          *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: orbital_index_radix, max_orbitals
      use debug_c
      use prnt_c
      use orb_c
      use buffer_c
      use eigv_c
      use ris_c
      use mpi_c,            only: mpi_rank, mpi_process_count
      use alcbuf_i
      use ispar_i
      use itjpo_i
      use convrt_i, only: convrt
      use rkco_gg_i, only: rkco_gg
      implicit none (type, external)
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      integer, intent(in) :: doit
      external :: cord
      integer, parameter :: key = orbital_index_radix
      character(len=11) :: cnum
      real(real64) :: evpair(nvec)
      real(real64) :: contri, contrik1, vcoeff, radial_k1, radial_sms,  &
         pair_scale, term_k1, term_sms
      integer :: j,k,loc,ic,ir,itjpoc,incor,lcnum, iia, iib, iic, iid
      integer :: i, occ_ic, occ_ir, ndiff

      incor = 1
      call alcbuf (1)

      do ic = mpi_rank + 1, csf_count, mpi_process_count
         if (mpi_rank == 0 .and. mod(ic,100) .eq. 1) then
            call convrt (ic,cnum,lcnum)
            print *, 'Column '//cnum(1:lcnum)//' complete on root;'
         endif

         itjpoc = itjpo (ic)
         do ir = ic,csf_count
            if (itjpo(ir) .eq. itjpoc) then
               if (ispar(ic) .ne. ispar(ir)) cycle
               ndiff = 0
               do i = 1,orbital_count
                  occ_ic = csf_subshell_occupations(i,ic)
                  occ_ir = csf_subshell_occupations(i,ir)
                  if (occ_ic .eq. occ_ir) cycle
                  if (iabs(occ_ic - occ_ir) .gt. 2) then
                     ndiff = 5
                     exit
                  endif
                  ndiff = ndiff + iabs(occ_ic - occ_ir)
                  if (ndiff .gt. 4) exit
               end do
               if (ndiff .gt. 4 .or. mod(ndiff,2) .ne. 0) cycle

               angular_coefficient_count = 0
               call rkco_gg (ic, ir, cord, incor, 1)
               if (angular_coefficient_count .eq. 0) cycle
               if (ir .eq. ic) then
                  pair_scale = 1.0_real64
               else
                  pair_scale = 2.0_real64
               endif
               do j = 1,nvec
                  loc = (j-1)*csf_count
                  evpair(j) = pair_scale * csf_mixing_coefficients(ic+loc) * csf_mixing_coefficients(ir+loc)
               end do
               do k = 1,angular_coefficient_count
                  vcoeff = angular_coefficients(k)
                  if (abs (vcoeff) .gt. cutoff) then
                     if (angular_coefficient_labels(5,k) .eq. 1) then
                        iia = angular_coefficient_labels(1,k)
                        iib = angular_coefficient_labels(2,k)
                        iic = angular_coefficient_labels(3,k)
                        iid = angular_coefficient_labels(4,k)
                        radial_k1 = vint(iia,iic) * vint(iib,iid)
                        radial_sms = (vint2(iia,iic) * vint(iib,iid) +   &
                           vint2(iib,iid) * vint(iia,iic))/2.0_real64
                        if (radial_k1 == 0.0_real64 .and.                     &
                            radial_sms == 0.0_real64) cycle
                        term_k1 = -vcoeff * radial_k1
                        term_sms = -vcoeff * radial_sms
                        if(doit.eq.1) write(51) ic,ir
                        if(doit.eq.1) then
                           write(51) vcoeff, ((iia*key + iic)*key+iib)*key+iid
                        endif
                        do j = 1,nvec
                           contrik1 = evpair(j) * term_k1
                           contri = evpair(j) * term_sms
                           smsc1(j) = smsc1(j) + contrik1
                           smsc2(j) = smsc2(j) + contri
                        end do
                     endif
                  endif
               end do
            endif
         end do
      end do

      call alcbuf (3)
      return
      end subroutine smsnew_mpi
