!***********************************************************************
!                                                                      *
      subroutine setham(jblock, myid, nprocs)
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:23:52   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: keyorb
      use memory_man
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      use hmat_c
      use mcpa_c
      use orb_c
      use pos_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i
      use rinti_i
      use fco_i
      use slater_i
      use gco_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: jblock
      integer , intent(in) :: myid
      integer , intent(in) :: nprocs
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
      character*6, parameter :: myname = 'SETHAM'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncfpat, ia, ir, itmp, idiag, km, k, k0, ib, nkjia, nkjib, kmin&
         , kmax, ndim, jblockt, ncft, ncoeff, ios, lab, ncontr, iclmndum, i, &
         loc, nfile, id, ic
      integer(long) :: lod
      real(double) :: diaa, coef, f0aa, fkaa, f0ab, gkab, tegral
      logical :: set
      character :: mcplab*3
      real(double), dimension(:), pointer :: coeff
      integer, dimension(:), pointer :: indx
!-----------------------------------------------
!
!     !*** needs EMT, IENDC.
!     POINTER (PNTEMT,EMT(1))
!     POINTER (PIENDC,IENDC(0:*))
!     POINTER (PNIROW,IROW(1))
!
!     !*** needs ncf, nw
!     POINTER (PNTRIQ,RIQDUMMY)
!
!
!     POINTER (pncfpast, ncfpast(1))
!     POINTER (pncminpast, ncminpast(1))
!     POINTER (pnevecpast, nevecpast(1))
!
!     POINTER (PCOEFF,COEFF(1))
!     !POINTER (PICLMN,ICLMNdum)
!     POINTER (PINDX,INDX(1))
!
!-----------------------------------------------------------------------
      ncfpat = ncfpast(jblock)
!=======================================================================
!   Accumulate diagonal terms that do not require MCP coefficients
!=======================================================================

!=======================================================================
!   Piece involving I(a,a) integrals
!=======================================================================

      do ia = 1, nw
         set = .false.
         do ir = myid + 1, ncf, nprocs
            itmp = iq(ia,ir + ncfpat)
            if (itmp <= 0) cycle
            !*** Occupation number not zero ...
            if (.not.set) then
               diaa = rinti(ia,ia,0)
               set = .true.
            endif
            ! IDIAG = IENDC(IR-1)+1
            idiag = iendc(ir)                    ! lower-triangle-by-rows mode
            emt(idiag) = emt(idiag) + itmp*diaa
         end do
      end do

!=======================================================================
!                    0
!   Piece involving F (a,a) integrals
!=======================================================================

      do ia = 1, nw
         set = .false.
         do ir = myid + 1, ncf, nprocs
            coef = fco(0,ir + ncfpat,ia,ia)
            if (coef == 0.d0) cycle
            !*** Angular coefficient not zero ...
            if (.not.set) then
               f0aa = slater(ia,ia,ia,ia,0)
               set = .true.
            endif
            ! IDIAG = IENDC(IR-1)+1
            idiag = iendc(ir)
            emt(idiag) = emt(idiag) + coef*f0aa
         end do
      end do

!=======================================================================
!                    k
!   Piece involving F (a,a) integrals
!=======================================================================

      km = 0
      k = 0
    6 continue
      k = k + 2
      do ia = 1, nw
         k0 = nkj(ia) - 1
         km = max0(k0,km)
         if (k > k0) cycle
         set = .false.

         do ir = myid + 1, ncf, nprocs
            coef = fco(k,ir + ncfpat,ia,ia)
            if (coef == 0.d0) cycle
            if (.not.set) then
               fkaa = slater(ia,ia,ia,ia,k)
               set = .true.
            endif
            ! IDIAG = IENDC(IR-1)+1
            idiag = iendc(ir)
            emt(idiag) = emt(idiag) + coef*fkaa
         end do
      end do
      if (k < km) go to 6

!=======================================================================
!                    0
!   Piece involving F (a,b) integrals
!=======================================================================

      do ia = 1, nw - 1
         do ib = ia + 1, nw
            set = .false.
            do ir = myid + 1, ncf, nprocs
               coef = fco(0,ir + ncfpat,ia,ib)
               if (coef == 0.d0) cycle
               if (.not.set) then
                  f0ab = slater(ia,ib,ia,ib,0)
                  set = .true.
               endif
            ! IDIAG = IENDC(IR-1)+1
               idiag = iendc(ir)
               emt(idiag) = emt(idiag) + coef*f0ab
            end do
         end do
      end do

!=======================================================================
!                    k
!   Piece involving G (a,b) integrals
!=======================================================================

      km = 0
      k = -1
   12 continue
      k = k + 1
      do ia = 1, nw - 1
         nkjia = nkj(ia)
         do ib = ia + 1, nw
            nkjib = nkj(ib)
            set = .false.
            if (nak(ia)*nak(ib) > 0) then
               kmin = abs((nkjia - nkjib)/2)
            else
               kmin = abs((nkjia - nkjib)/2) + 1
            endif
            if (mod(k - kmin,2) /= 0) cycle

            kmax = (nkjia + nkjib)/2
            km = max0(kmax,km)
            if (k<kmin .or. k>kmax) cycle

            do ir = myid + 1, ncf, nprocs
               coef = gco(k,ir + ncfpat,ia,ib)
               if (coef == 0.d0) cycle
               if (.not.set) then
                  gkab = slater(ia,ib,ib,ia,k)
                  set = .true.
               endif
               ! IDIAG = IENDC(IR-1)+1
               idiag = iendc(ir)
               emt(idiag) = emt(idiag) + coef*gkab
            end do
         end do
      end do
      if (k < km) go to 12

!=======================================================================
!   Local storage for reading mcpXXX files
!=======================================================================

      ndim = 1
      call alloc (coeff, ndim, 'COEFF', 'SETHAM' )
      call alloc (indx, ndim, 'INDX', 'SETHAM')

!=======================================================================
!   Accumulate one-electron terms that require MCP coefficients
!=======================================================================

      read (31) mcplab, jblockt, ncft, ncoeff

      if (jblockt /= jblock) then
         write (istde, *) myname, ': blk1=', jblockt, ' blk2=', jblock
         error stop
      endif
      if (ncft /= ncf) then
         write (istde, *) myname, ': ncf1 = ', ncft, ' ncf2 = ', ncf
         error stop
      endif

!=======================================================================
!   Loop over non-zero labels which have non-zero elements
!=======================================================================

      read (31, iostat=ios) lab, ncontr
      if (ios /= 0) error stop 'IOS .NE. 0 when reading LAB, NCONTR'
      do while(lab/=0 .or. ncontr/=0)

         !*** decode the label of I(ab)
         ia = mod(lab,key)
         ib = lab/key

         !*** Compute radial integral I(ab)
         tegral = rinti(ia,ib,0)

         ! Read column index, sparse matrix index, and coefficient
         ! for all contributions from this integral.
         if (ncontr > ndim) then
            call dalloc (coeff, 'COEFF', 'SETHAM')
            call dalloc (indx, 'INDX', 'SETHAM')
            ndim = ncontr
            call alloc (coeff, ndim, 'COEFF', 'SETHAM')
            call alloc (indx, ndim, 'INDX', 'SETHAM')
         endif
         read (31) (iclmndum,indx(i),coeff(i),i=1,ncontr)

         !*** Store all the contributions from this integral
         do i = 1, ncontr
            loc = indx(i)
            lod = loc
            if (lod > (nelmnt)) then
               write (6, *) '  Error in computing 1-e contribution'
               write (6, *) '  LOC = ', loc, '  NELMNT = ', (nelmnt)
               error stop
            endif
            emt(loc) = emt(loc) + tegral*coeff(i)
         end do

         read (31, iostat=ios) lab, ncontr
         if (ios == 0) cycle
         error stop 'IOS .NE. 0 when reading LAB, NCONTR'
      end do

!=======================================================================
!   Accumulate two-electron terms that require MCP coefficients
!=======================================================================

      do nfile = 32, 32 + kmaxf
         k = nfile - 32

         read (nfile) mcplab, jblockt, ncft, ncoeff

         if (jblockt /= jblock) then
            write (istde, *) myname, ': blk3=', jblockt, ' blk4=', jblock
            error stop
         endif
         if (ncft /= ncf) then
            write (istde, *) myname, ': ncf3 = ', ncft, ' ncf4 = ', ncf
            error stop
         endif

!=======================================================================
!   Loop over non-zero labels which have non-zero elements
!=======================================================================

         read (nfile, iostat=ios) lab, ncontr
         if (ios /= 0) error stop 'IOS .NE. 0 when reading LAB, NCONTR 2'
         do while(lab/=0 .or. ncontr/=0)

            !                         k
            !*** decode the label of R (abcd)
            id = mod(lab,key)
            lab = lab/key
            ib = mod(lab,key)
            lab = lab/key
            ic = mod(lab,key)
            ia = lab/key

            !*** Compute radial integral
            tegral = slater(ia,ib,ic,id,k)

            ! Read column index, sparse matrix index, and coefficient
            ! for all contributions from this integral.
            if (ncontr > ndim) then
               call dalloc (coeff, 'COEFF', 'SETHAM')
               call dalloc (indx, 'INDX', 'SETHAM')
               ndim = ncontr
               call alloc (coeff, ndim, 'COEFF', 'SETHAM')
               call alloc (indx, ndim, 'INDX', 'SETHAM')
            endif
            read (nfile) (iclmndum,indx(i),coeff(i),i=1,ncontr)

            !*** Store all the contributions from this integral
            do i = 1, ncontr
               loc = indx(i)
               lod = loc
               if (lod > nelmnt) then
                  write (6, *) '  Error in computing 2-e contribution'
                  write (6, *) '  LOC = ', loc, '  NELMNT = ', nelmnt
                  error stop
               endif
               emt(loc) = emt(loc) + tegral*coeff(i)
            end do

            read (nfile, iostat=ios) lab, ncontr
            if (ios == 0) cycle
            error stop 'IOS .NE. 0 when reading LAB, NCONTR 2'
         end do
      end do

!=======================================================================
!   Deallocate local storage
!=======================================================================

      call dalloc (coeff, 'COEFF', 'SETHAM')
      call dalloc (indx, 'INDX', 'SETHAM')

      return
      end subroutine setham
