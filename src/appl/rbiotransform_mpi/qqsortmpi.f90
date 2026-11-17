!***********************************************************************
!                                                                      *
      subroutine qqsort(nfile, number, kstart, name, kamax)
!                                                                      *
!     The list of unique integrals (j,i) is formed in the order of     *
!     increasing symmetry, i.e. with j .le. i.                         *
!     With each integral (INTGRL) there is a pointer (INTPTR)          *
!     that points to the last element in the array of coefficients     *
!     for that integral.                                               *
!                                                                      *
!     On Exit                                                          *
!                                                                      *
!     INTGRL -- array of integrals (packed form of j,i)                *
!     NINTG  -- number of integrals                                    *
!     INTPTR -- array of pointers to last element in list of coeff.    *
!     CNN    -- array of coefficients                                  *
!     NCOEFF -- number of coefficients                                 *
!     JANN   -- row of matrix                                          *
!     JBNN   -- column of matrix                                       *
!                                                                      *
!     Written by Per Jonsson                                           *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use memory_man, only: alloc, dalloc, double
      use orb_c, only: csf_count, double, orbital_count
      use default_c, only: double, inpci, ndef, ndump
      use mcpdata_c, only: cnn, double, intgrl, intptr, jann, jbnn, ncoeff, nintg
      use mpi_c, only: double, mpi_process_count, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nfile
      integer, intent(in) :: number
      integer  :: kstart
      integer, intent(in) :: kamax
      character  :: name*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
      integer, parameter :: nf = 200
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, l, ir, ja, jb, int, j, intt
      real(real64) :: cn
!-----------------------------------------------
!
      rewind (nfile + 80)

      ncoeff = number
!
!   Sort the list
!
!   Allocate storage for all required arrays these arrays are then deallocated
!   in mcp
!
      call alloc (jann,   ncoeff, "JANN",   "QQSORT")
      call alloc (jbnn,   ncoeff, "JBNN",   "QQSORT")
      call alloc (intgrl, ncoeff, "INTGRL", "QQSORT")
      call alloc (cnn,    ncoeff, "CNN",    "QQSORT")
      call alloc (intptr, ncoeff, "INTPTR", "QQSORT")
!
!   Read arrays into memory from NFILE
!
      do i = 1, ncoeff
         read (nfile + 80) jann(i), jbnn(i), intgrl(i), cnn(i)
      end do
!
!   Sort INTGRL into ascending order using the heapsort algorithm;
!   (Numerical recepies page 231.) move the associated members of
!   of JANN and JBNN in the same
!   manner;
!
      if (ncoeff > 1) then
!
         l = ncoeff/2 + 1
         ir = ncoeff
    2    continue
         if (l > 1) then
            l = l - 1
            ja = jann(l)
            jb = jbnn(l)
            int = intgrl(l)
            cn = cnn(l)
         else
            ja = jann(ir)
            jb = jbnn(ir)
            int = intgrl(ir)
            cn = cnn(ir)
            jann(ir) = jann(1)
            jbnn(ir) = jbnn(1)
            intgrl(ir) = intgrl(1)
            cnn(ir) = cnn(1)
            ir = ir - 1
            if (ir == 1) then
               jann(1) = ja
               jbnn(1) = jb
               intgrl(1) = int
               cnn(1) = cn
               goto 4
            end if
         end if
         i = l
         j = l + l
    3    continue
         if (j <= ir) then
            if (j < ir) then
               if (intgrl(j) < intgrl(j+1)) j = j + 1
            end if
            if (int < intgrl(j)) then
               jann(i) = jann(j)
               jbnn(i) = jbnn(j)
               intgrl(i) = intgrl(j)
               cnn(i) = cnn(j)
               i = j
               j = j + j
            else
               j = ir + 1
            end if
            goto 3
         end if
         jann(i) = ja
         jbnn(i) = jb
         intgrl(i) = int
         cnn(i) = cn
         goto 2
      end if
!

!   Sorting complete; close the file
!
!ww    4 CLOSE (80+NFILE)
    4 continue
      nintg = 1
      intt = intgrl(1)
!
      do i = 1, ncoeff
         if (intgrl(i) == intt) cycle
         intptr(nintg) = i - 1
         nintg = nintg + 1
         intt = intgrl(i)
      end do

      intptr(nintg) = ncoeff
!
      do i = 1, nintg
         intgrl(i) = intgrl(intptr(i))
      end do
!
!  If output option is set dump the data on file
!
      if (ndump == 1) then
!
!  If first set of data open the file and print
!  some data to later be able to identify the file
!
!        IF (NFILE.EQ.1.AND.IBLK.EQ.1) THEN
!         print *, ' open sorted ang. file .TB(NF)', NF
!          J = INDEX(NAME,' ')
!          OPEN (UNIT = NF,FILE=NAME(1:J-1)//'.TB',
!     :          STATUS='UNKNOWN',FORM='UNFORMATTED')
!          REWIND (NF)
!        ENDIF
         if (nfile == 1) write (nf) csf_count,orbital_count,kamax,mpi_process_count,mpi_rank
!
!  Print out angular data for this kappa
!
         write (nf) nintg, ncoeff
         do i = 1, nintg
            write (nf) intgrl(i), intptr(i)
         end do
         do i = 1, ncoeff
            write (nf) cnn(i), jann(i), jbnn(i)
         end do
      end if
      call dalloc (jann,   "JANN",   "QQSORT")
      call dalloc (jbnn,   "JBNN",   "QQSORT")
      call dalloc (intgrl, "INTGRL", "QQSORT")
      call dalloc (cnn,    "CNN",    "QQSORT")
      call dalloc (intptr, "INTPTR", "QQSORT")
!
!  Has all data been processed? If so close
!  the file
!
!     IF (NFILE.EQ.KAMAX) CLOSE (NF)

!  Debug output
!
!      IF (KSTART.EQ.1) THEN
!        NNNN = 120
!      ELSE
!        NNNN = 140
!      ENDIF
!      WRITE(NNNN+NFILE,*) ' Data for '
!      WRITE(NNNN+NFILE,*) NINTG, ' Integrals'
!      DO I = 1,NINTG
!        IA = INTGRL(I)/KEY
!        IB = MOD(INTGRL(I),KEY)
!        WRITE(NNNN+NFILE,302) NP(IA),NH(IA),NP(IB),NH(IB),INTPTR(I)
!      ENDDO
!      WRITE(NNNN+NFILE,*) NCOEFF,' Coefficients'
!      WRITE(NNNN+NFILE,'(F12.8,2I6)')
!     :     (CNN(I),JANN(I),JBNN(I),I=1,NCOEFF)
!
      return
!
  301 format(" T_[",1i2,",",1i2,"]"," (",1i2,1a2,",",1i2,1a2,") = ",1p,d19.12)
  302 format(" (",1i2,1a2,",",1i2,1a2,")",i6)
      return
!
      end subroutine qqsort
