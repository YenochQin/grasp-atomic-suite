!***********************************************************************
!                                                                      *
      subroutine mctout(iopar, jkp, name)
!                                                                      *
!   This routine loads  coefficients with parity and  rank specified   *
!   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
!   (+/- 1) and is determined from the sign of  KP(JKP).               *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:29:28   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: keyorb, nnnw
      use memory_man
      use blk_c
      use debug_c,         only: ldbpa
      use decide_c
      use foparm_c
      use orb_c
      use offd_c
      use osc_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use angdata_i
      use itjpo_i
      use oneparticlejj_i
      use trsort_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: iopar
      integer  :: jkp
      character  :: name(2)*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: cutoff = 1.0d-10
      integer, parameter :: nfile = 93
      integer, parameter :: nfile1 = 237
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nfile2, iblki, iblkf, nmct, nlabel, ncfi0, ncff0, ic, ir, ncr&
         , ia, ib, newsiz, i
      logical :: avail
      integer, dimension(:), pointer :: label
      real(double), dimension(:), pointer:: coeff
      real(double), dimension(nnnw) :: tshell
!-----------------------------------------------
!
!
!
! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
!
!   Check if angular data is available on file
!
      nfile2 = nfile1 + jkp
      call angdata (name, avail, jkp, nfile2)

!
!   If angular data is not available open the scratch file to store the
!   coefficients; position file
!   to beginning
!
!
      lk = abs(kp(jkp))
      iopar = isign(1,kp(jkp))
      if (avail) return
      write (6, *) 'LK,IOPAR,from MCTOUT'
      write (6, *) lk, iopar
!
! Start of the block loops
      do iblki = 1, nblocki
         do iblkf = 1, nblockf
!        OPEN (NFILE,STATUS = 'new', FORM = 'UNFORMATTED')
!        Sometimes, when there has been an error, status need to be "unknown"
            open(nfile,status='unknown',form='UNFORMATTED')
            nmct = 0
!
!
!   If angular data is not available
!   Generate MCT coefficients for given rank/parity combination and
!   store them by CSF on NFILE
!
!
!   Allocate storage to buffer arrays
!
            nlabel = 32
            call alloc (label, nlabel, 'LABEL', 'MCTOUT')
            call alloc (coeff, nlabel, 'COEFF', 'MCTOUT')

            if (iblki == 1) then
               ncfi0 = 1
            else
               ncfi0 = ncfi(iblki-1) + 1
            endif
!
            if (iblkf == 1) then
               ncff0 = ncfi(nblocki) + 1
            else
               ncff0 = ncfi(nblocki) + ncff(iblkf-1) + 1
            endif
            do ic = ncfi0, ncfi(iblki)
               do ir = ncff0, ncfi(nblocki) + ncff(iblkf)
!
!           IR = IC
!
                  ncr = 0

!
!   In many case one is interested only in M1 and E2 transitions between
!   levels with different J values. If this is the case then the do check
!   on the J quantum numbers of the CSFs before calling TNSRJJ.
!
                  if (kp(jkp)==1 .and. noffd1==1) then
                     if (itjpo(ic) == itjpo(ir)) cycle
                  endif
                  if (kp(jkp)==2 .and. noffd2==1) then
                     if (itjpo(ic) == itjpo(ir)) cycle
                  endif
!          if(ispar(ic)*ispar(ir)*iopar.ne.1.
!    &        or.itrig(itjpo(ic),itjpo(ir),2*lk+1).ne.1) go to 13
!          if(ichkq1(IC,IR).eq.0) go to 13
                  call oneparticlejj(lk,iopar,ic,ir,ia,ib,tshell)
                  if (ia /= 0) then
                     if (ia == ib) then
                        do ia = 1, nw
                           if (abs(tshell(ia)) <= cutoff) cycle
                           ncr = ncr + 1
                           if (ncr > nlabel) then
                              newsiz = 2*nlabel
                              call ralloc (label,  newsiz, 'LABEL', 'MCTOUT')
                              call ralloc (coeff,  newsiz, 'COEFF', 'MCTOUT')
                              nlabel = newsiz
                           endif
                           label(ncr) = ia*keyorb + ia
                           coeff(ncr) = tshell(ia)
                        end do
                     else
                        if (abs(tshell(1)) > cutoff) then
                           ncr = ncr + 1
                           if (ncr > nlabel) then
                              newsiz = 2*nlabel
                              call ralloc (label,  newsiz, 'LABEL', 'MCTOUT')
                              call ralloc (coeff,  newsiz, 'COEFF', 'MCTOUT')
                              nlabel = newsiz
                           endif
                           label(ncr) = ia*keyorb + ib
                           coeff(ncr) = tshell(1)
                        endif
                     endif
                  endif
                  if (ncr <= 0) cycle
                  write (nfile) ic - ncfi0 + 1, ir - ncff0 + 1, ncr
!             WRITE (NFILE) IC,IR,NCR
                  write (nfile) (label(i),coeff(i),i=1,ncr)
                  nmct = nmct + ncr

!
               end do
            end do
!
!   Deallocate storage for buffer arrays
!
            call dalloc (label,  'LABEL', 'MCTOUT')
            call dalloc (coeff,  'COEFF', 'MCTOUT')
!
            write (*, 301) nmct, lk, iopar
!
!   Sort the MCT coefficients by integral labels
!
            call trsort (name, nfile, nfile2, ldbpa(2), jkp, iblki, iblkf)
            close(nfile, status='delete')
! end of the loops for blocks
         end do
      end do
!
!   Read the data back as required by OSCL conventions
!
      rewind (nfile2)
      return
!
  301 format(/,/,/,1x,i8,' MCT coefficients generated for rank ',i2,&
         ' and parity ',i2,/,/)
      return
!
      end subroutine mctout
