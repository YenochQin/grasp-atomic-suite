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
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use orb_c
      use buffer_c
      use eigv_c
      use ris_c
      use mpi_c,            only: mpi_rank, mpi_process_count
      use alcbuf_i, only: alcbuf
      use ispar_i, only: ispar
      use itjpo_i, only: itjpo
      use convrt_i, only: convrt
      use rkco_gg_i, only: rkco_gg
      use cord_i, only: cord
      implicit none (type, external)
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      integer, intent(in) :: doit

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
         if (mpi_rank == 0 .and. mod(ic,100) == 1) then
            call convrt (ic,cnum,lcnum)
            print *, "Column "//cnum(1:lcnum)//" complete on root;"
         end if

         itjpoc = itjpo (ic)
         do ir = ic,csf_count
            if (itjpo(ir) == itjpoc) then
               if (ispar(ic) /= ispar(ir)) cycle
               ndiff = 0
               do i = 1,orbital_count
                  occ_ic = csf_subshell_occupations(i,ic)
                  occ_ir = csf_subshell_occupations(i,ir)
                  if (occ_ic == occ_ir) cycle
                  if (iabs(occ_ic - occ_ir) > 2) then
                     ndiff = 5
                     exit
                  end if
                  ndiff = ndiff + iabs(occ_ic - occ_ir)
                  if (ndiff > 4) exit
               end do
               if (ndiff > 4 .or. mod(ndiff,2) /= 0) cycle

               angular_coefficient_count = 0
               call rkco_gg (ic, ir, cord, incor, 1)
               if (angular_coefficient_count == 0) cycle
               if (ir == ic) then
                  pair_scale = 1.0_real64
               else
                  pair_scale = 2.0_real64
               end if
               do j = 1,nvec
                  loc = (j-1)*csf_count
                  evpair(j) = pair_scale * csf_mixing_coefficients(ic+loc) * csf_mixing_coefficients(ir+loc)
               end do
               do k = 1,angular_coefficient_count
                  vcoeff = angular_coefficients(k)
                  if (abs (vcoeff) > cutoff) then
                     if (angular_coefficient_labels(5,k) == 1) then
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
                        if(doit==1) write(51) ic,ir
                        if(doit==1) then
                           write(51) vcoeff, ((iia*key + iic)*key+iib)*key+iid
                        end if
                        do j = 1,nvec
                           contrik1 = evpair(j) * term_k1
                           contri = evpair(j) * term_sms
                           smsc1(j) = smsc1(j) + contrik1
                           smsc2(j) = smsc2(j) + contri
                        end do
                     end if
                  end if
               end do
            end if
         end do
      end do

      call alcbuf (3)
      return
      end subroutine smsnew_mpi
