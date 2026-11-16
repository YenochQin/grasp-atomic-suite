!***********************************************************************
!                                                                      *
      subroutine trsort(name, nfile, nfile2, lprint, jkp, iblki, iblkf)
!                                                                      *
!   Routine to sort angular coefficients into list based on integral   *
!   labels rather than CSF.  A tree sort is used. To save space, the   *
!   CSF pair labels and coefficients  are not read into arrays until   *
!   the sorting has been done.                                         *
!                                                                      *
!   Call(s) to: [OSCL92]: ALCLLA, ALCNMA.                              *
!                                                                      *
!   Original author(s) unknown. Modifications for dynamic storage by   *
!   Farid A. Parpia.                                                   *
!                                                                      *
!                                           Last update: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use memory_man
      use default_c, only: double, inpci, ndef, ndump
      use orb_c
      use osc_c, only: nkp, kp
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alclla_i, only: alclla
      use alcnma_i, only: alcnma
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nfile
      integer, intent(in) :: nfile2
      integer  :: jkp
      integer, intent(in) :: iblki
      integer, intent(in) :: iblkf
      logical, intent(in) :: lprint
      character, intent(in) :: name(2)*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nca = 65536
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, pointer :: ilab(:), iright(:), ileft(:), ibeg(:), iptr(:), &
                    ilast(:), iblint(:), iptcsf(:), lblint(:), jlabl(:), isldr(:), isldr1(:)
      real(real64), pointer :: xsldr(:), xl(:)
      integer :: nmcp, nint, nlabel, lldim, nmdim, ir, is, ni, i, m, j, &
                 icount, jlab, k, l, imcp, j1, j2, mlr, mup, ints, ia, ib
      real(real64) :: x
      logical :: first
      character(len=2) :: s(-9:9)
!-----------------------------------------------
!
!
!     POINTER (PIBEG,IBEG(1)),(PILAB,ILAB(1)),(PILAST,ILAST(1)),
!    :        (PILEFT,ILEFT(1)),(PIPTCS,IPTCSF(1)),(PIRIGH,IRIGHT(1)),
!    :        (PLBLIN,LBLINT(1))
!     POINTER (PIPTR,IPTR(1))
!
!
      s((-9)) = "-9"
      s((-8)) = "-8"
      s((-7)) = "-7"
      s((-6)) = "-6"
      s((-5)) = "-5"
      s((-4)) = "-4"
      s((-3)) = "-3"
      s((-2)) = "-2"
      s((-1)) = "-1"
      s(0) = "+0"
      s(1) = "+1"
      s(2) = "+2"
      s(3) = "+3"
      s(4) = "+4"
      s(5) = "+5"
      s(6) = "+6"
      s(7) = "+7"
      s(8) = "+8"
      s(9) = "+9"
!
!   Position file at beginning of list of integrals
!
      rewind (nfile)
!
!   Initialize
!
      first = .true.
      nmcp = 0
      nint = 0
!
!   Initial allocation of storage to local arrays
!
      nlabel = 1
      call alloc (jlabl, nlabel, "JLABL", "TRSORT")
      call alloc (xl, nlabel, "XL", "TRSORT")
!
      call alclla (ibeg, ilab, ilast, ileft, iptcsf, iright, lblint, &
          lldim, 1)
      call alcnma (iptr, isldr, isldr1, xsldr, nmdim, 1)
!
!   Now the rest of the elements
!
    1 continue
      read (nfile, end=12) ir, is, ni
      if (ni > nlabel) then
         call ralloc (jlabl, ni, "JLABL", "TRSORT")
         call ralloc (xl, ni, "XL", "TRSORT")
         nlabel = ni
      end if
      read (nfile, end=99, err=99) (jlabl(i),xl(i),i=1,ni)
      if (ir==0 .or. is==0 .or. ni==0) goto 1
!
      if (first) then
!
!   List is empty
!
         m = 0
         j = 0
!
!   Set up list pointers and insert first element
!
         icount = 0
    3    continue
         icount = icount + 1
         if (icount > ni) goto 1
         if (jlabl(icount) == 0) goto 3
         ilab(1) = jlabl(icount)
!
         iright(1) = 0
         ileft(1) = 0
         ibeg(1) = 1
         iptr(1) = 0
         ilast(1) = 0
!
         m = 1
         j = 1
!
         first = .false.
!
      else
!
         icount = 0
!
      end if
!
!   Sort integral list using tree sort
!
!   Take next nonzero element
!
    4 continue
      icount = icount + 1
      if (icount > ni) goto 1
      jlab = jlabl(icount)
      if (jlab == 0) goto 4
      x = xl(icount)
!
      m = m + 1
      if (m > nmdim) call alcnma (iptr, isldr, isldr1, xsldr, nmdim, 2)
      i = 1
!
!   Search for place in tree
!
    5 continue
      if (jlab - ilab(i) > 0) goto 8
      if (jlab - ilab(i) == 0) goto 10
      k = iright(i)
      if (k /= 0) goto 7
      j = j + 1
      if (j > lldim) call alclla (ibeg, ilab, ilast, ileft, iptcsf, &
                    iright, lblint, lldim, 2)
      iright(i) = j
      goto 9
    7 continue
      i = k
      goto 5
    8 continue
      k = ileft(i)
      if (k /= 0) goto 7
      j = j + 1
      if (j > lldim) call alclla (ibeg, ilab, ilast, ileft, iptcsf, &
                          iright, lblint, lldim, 2)
      ileft(i) = j
!
!   When found, update list.
!
    9 continue
      ilast(j) = i
      iright(j) = 0
      ileft(j) = 0
      ibeg(j) = m
      ilab(j) = jlab
      iptr(m) = 0
      goto 4
   10 continue
      k = ibeg(i)
      l = k
      k = iptr(l)
      do while(k /= 0)
         l = k
         k = iptr(l)
      end do
      iptr(l) = m
      iptr(m) = 0
      goto 4
!
!   The end of the CSF-based file has been reached
!
   12 continue
      if (.not.(first .or. m==0)) then
!
!  Sort is complete. Unpack list
!
         nmcp = m
         nint = j
         l = 0
         m = 0
         i = 1
!
!   Search for smallest element
!
   13    continue
         k = iright(i)
         do while(k /= 0)
            i = k
            k = iright(i)
         end do
!
!   Insert in sorted list
!
   14    continue
         if (ilab(i) == 0) goto 16
         l = l + 1
         lblint(l) = ilab(i)
         k = ibeg(i)
!
!   Copy list of pointers to CSF/coefficients into new list
!
         m = m + 1
         isldr(m) = k
         k = iptr(k)
         do while(k /= 0)
            m = m + 1
            isldr(m) = k
            k = iptr(k)
         end do
         iptcsf(l) = m
         ilab(i) = 0
!
!   Next smallest element is on left of last element
!
         k = ileft(i)
         if (k == 0) goto 16
         i = k
         goto 13
!
!   If no element on left, next smallest is previous element
!
   16    continue
         i = ilast(i)
         if (i /= 0) goto 14
!
!   List is unpacked. Invert CSF/coefficient pointer list to give
!   position list for CSF/coefficients as they are read in
!
         do i = 1, nmcp
            k = isldr(i)
            iptr(k) = i
         end do
!
!   Now read CSF pairs and coefficients into correct positions in
!   sorted list
!
         imcp = 0
         rewind (nfile)
   18    continue
         read (nfile, end=20) ir, is, ni
         if (ir==0 .or. is==0 .or. ni==0) goto 18
         read (nfile) (jlabl(i),xl(i),i=1,ni)
         do i = 1, ni
            if (jlabl(i) == 0) cycle
            imcp = imcp + 1
            k = iptr(imcp)
            isldr(k) = ir
            isldr1(k) = is
!           ISLDR(K) = IR*NCA+IS
            xsldr(k) = xl(i)
         end do
         goto 18
!
!   The integral-based list is completely known
!
      end if
   20 continue
      rewind (nfile)
!
!  If first set of data open the file and print
!  some data to later be able to identify the file
!
      if (iblki==1 .and. iblkf==1) then
         j1 = index(name(1)," ")
         j2 = index(name(2)," ")
         open(unit=nfile2,file=name(1)(1:j1-1)//"."//name(2)(1:j2-1)//"."//s(&
            kp(jkp))//"T", status="UNKNOWN", form="UNFORMATTED", position=&
            "asis")
      end if
      write (nfile2) iblki,iblkf,orbital_count,nkp,mpi_process_count,mpi_rank
      write (nfile2) nint
      if (nmcp /= 0) then
         if (lprint) write (99, 301)
         mlr = 1
         do i = 1, nint
            mup = iptcsf(i)
            ints = mup - mlr + 1
            write (nfile2) lblint(i), ints
            if (lprint) then
               ia = lblint(i)/key
               ib = mod(lblint(i),key)
               write (99, 302) orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                  & orbital_symmetry_labels(ib)
            end if
            write (nfile2) (isldr(m),isldr1(m),xsldr(m),m=mlr,mup)
            if (lprint) then
               do m = mlr, mup
!                  IS = MOD (ISLDR(M),NCA)
                  is = isldr1(m)
!                 IR = ISLDR(M)/NCA
                  ir = isldr(m)
                  write (99, 303) ir, is, xsldr(m)
               end do
            end if
            mlr = mup + 1
         end do
!
      end if
!
!   Deallocate storage for local arrays
!
      call dalloc (jlabl, "JLABJ", "TRSORT")
      call dalloc (xl, "XL", "TRSORT")
      call alclla (ibeg, ilab, ilast, ileft, iptcsf, iright, lblint, &
                    lldim, 3)
      call alcnma (iptr, isldr, isldr1, xsldr, nmdim, 3)
!
      return
!
!   Error handling
!
   99 continue
      write (6, *) "TRSORT: Error reading CSF-based file."
      stop
!
  301 format(/,/,"  k"/," d  (rs)  Coefficients:"/,"  ab"/,/,/,&
         "    a     b          r         s  Coefficient"/)
  302 format(2(2x,i2,a2))
  303 format(14x,1i6,2x,1i6,2x,1p,1d22.15)
      return
!
      end subroutine trsort
