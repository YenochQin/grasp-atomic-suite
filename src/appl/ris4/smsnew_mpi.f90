!***********************************************************************
!                                                                      *
      subroutine smsnew_mpi(doit,vint,vint2)
!                                                                      *
!   MPI version of SMSNEW. Each rank handles a strided subset of the   *
!   outer IC loop and can optionally write a per-rank angular-data      *
!   fragment that the root process later merges into name.ITB.          *
!                                                                      *
!***********************************************************************
      use vast_kind_param,  only: double
      use parameter_def,    only: keyorb, nnnw
      use debug_c
      use prnt_c
      use orb_c
      use buffer_c
      use eigv_c
      use ris_c
      use mpi_c,            only: myid, nprocs
      use alcbuf_i
      use ispar_i
      use itjpo_i
      use convrt_i, only: convrt
      use rkco_gg_i, only: rkco_gg
      implicit none
      real(double), dimension(nnnw,nnnw), intent(in) :: vint, vint2
      integer, intent(in) :: doit
      external :: cord
      integer, parameter :: key = keyorb
      character*11 :: cnum
      real(double), dimension(nvec) :: evpair
      real(double) :: contri, contrik1, vcoeff, radial_k1, radial_sms,  &
         pair_scale, term_k1, term_sms
      integer :: j,k,loc,ic,ir,itjpoc,incor,lcnum, iia, iib, iic, iid
      integer :: i, occ_ic, occ_ir, ndiff

      incor = 1
      call alcbuf (1)

      do ic = myid + 1, ncf, nprocs
         if (myid == 0 .and. mod(ic,100) .eq. 1) then
            call convrt (ic,cnum,lcnum)
            print *, 'Column '//cnum(1:lcnum)//' complete on root;'
         endif

         itjpoc = itjpo (ic)
         do ir = ic,ncf
            if (itjpo(ir) .eq. itjpoc) then
               if (ispar(ic) .ne. ispar(ir)) cycle
               ndiff = 0
               do i = 1,nw
                  occ_ic = iqa(i,ic)
                  occ_ir = iqa(i,ir)
                  if (occ_ic .eq. occ_ir) cycle
                  if (iabs(occ_ic - occ_ir) .gt. 2) then
                     ndiff = 5
                     exit
                  endif
                  ndiff = ndiff + iabs(occ_ic - occ_ir)
                  if (ndiff .gt. 4) exit
               end do
               if (ndiff .gt. 4 .or. mod(ndiff,2) .ne. 0) cycle

               nvcoef = 0
               call rkco_gg (ic, ir, cord, incor, 1)
               if (nvcoef .eq. 0) cycle
               if (ir .eq. ic) then
                  pair_scale = 1.0d00
               else
                  pair_scale = 2.0d00
               endif
               do j = 1,nvec
                  loc = (j-1)*ncf
                  evpair(j) = pair_scale * evec(ic+loc) * evec(ir+loc)
               end do
               do k = 1,nvcoef
                  vcoeff = coeff(k)
                  if (abs (vcoeff) .gt. cutoff) then
                     if (label(5,k) .eq. 1) then
                        iia = label(1,k)
                        iib = label(2,k)
                        iic = label(3,k)
                        iid = label(4,k)
                        radial_k1 = vint(iia,iic) * vint(iib,iid)
                        radial_sms = (vint2(iia,iic) * vint(iib,iid) +   &
                           vint2(iib,iid) * vint(iia,iic))/2.0d00
                        if (radial_k1 == 0.0d00 .and.                     &
                            radial_sms == 0.0d00) cycle
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
