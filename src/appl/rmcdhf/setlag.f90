!***********************************************************************
!                                                                      *
      subroutine setlag(eol)
!                                                                      *
!   Sets up the data structure  pertaining to the Lagrange multipli-   *
!   ers  on the first entry;  on subsequent calls it  determines new   *
!   estimates for the multipliers.                                     *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, QUAD, RINTI,                           *
!               [RSCF92]: DACON, SETCOF, XPOT,  YPOT.                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 18 Dec 1992   *
!   MPI version by Xinghong He              Last update: 03 Aug 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:25:08   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: keyorb
      use memory_man
      use orba_c
      use core_c
      use def_c
      use fixd_c
      use grid_c, only: n, rpor
      use lagr_c
      use orb_c
      use mpi_s
      use pote_c, only: yp, xp, xq
      use scf_c
      use tatb_c, only: mtp, ta
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcof_i
      use ypot_i
      use xpot_i
      use dacon_i
      use quad_i
      use rinti_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical  :: eol
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: p001 = 1.0d-01
      integer, parameter :: key = keyorb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: itwice, liraw, li, lip1, naklj, ljraw, lj, ieccli, l1, l2, &
         jlast, mlast, m, j, i
      real(double), dimension(nnnp) :: ypj, ypm, xpj, xpm, xqj, xqm
      real(double) :: eps, ucfj, ucfm, result, rijm, qdif, obqdif, obqsum
      logical :: first, fixli, fixlj, fulli, fullj, vli, vlj
!-----------------------------------------------
!
      data first/ .true./
!
!-----------------------------------------------------------------------

      if (first) then

!=======================================================================
!   Determine the total number of Lagrange multipliers and store
!   their indeces in IECC(1:NEC). Memories are allocated for IECC
!   and ECV. The only outputs are NEC and IECC(). This is done only
!   for the first time this routine is called.
!
!   Implementation changes:
!   The "DO LIraw ..." is executed twice (if NEC obtained in the
!   first time > 0) so the memory allocation is done only once.
!
!   This part is not distributed.
!=======================================================================

         eps = accy*0.01d0      ! criterion to see if an orb is occupied
         do itwice = 1, 2
            nec = 0
         do ljraw = ncore+1, nw
               lj = iorder(ljraw)
            naklj = nak(lj)
            vlj = .not. lfix(lj)
               fullj = abs ( ucf(lj)-dble (nkj(lj)+1) ) .lt. eps
            do liraw = 1, ljraw-1
               li = iorder(liraw)
               vli = .not. lfix(li)
               fulli = abs ( ucf(li)-dble (nkj(li)+1) ) .lt. eps
               if  (nak(li) .eq. naklj) then
                  if  (vli  .or. vlj ) then       ! at least one orbital varied
!                                                   but not both (varied and full)
                     if (.not. ((vli .and. vlj) .and. (fulli .and. fullj))) then 
                  nec = nec + 1
                  !*** Encode index at 2nd round ***
                  if (itwice == 2) iecc(nec) = li + key * lj
               endif
                  endif
               endif
            enddo
         enddo

            if (itwice==1 .and. nec>0) then
               call alloc (ecv, nec, 'ECV', 'SETLAG')
               call alloc (iecc, nec, 'IECC', 'SETLAG')
            else
               exit
            endif
         end do                                  !itwice

!=======================================================================
!   Print information about Lagrange multipliers
!=======================================================================

         if (myid == 0) then
            if (nec == 0) then
               write (*, 302)
            else
               write (*, 304)
               do li = 1, nec
               !*** Decode index ***
                  ieccli = iecc(li)
                  l1 = ieccli/key
                  l2 = ieccli - key*l1
                  write (*, 305) np(l2), nh(l2), np(l1), nh(l1)
               end do
            endif
         endif
         first = .false.
      endif

!     Lagrange multipliers need to be computed also on the first call
!     RETURN

      if (nec == 0) return
      if (myid == 0) write (*, 306)
      jlast = 0
      mlast = 0

      do li = 1, nec
         !*** Decode index ***
         ieccli = iecc(li)
         m = ieccli/key
         j = ieccli - key*m
!
         if (j /= jlast) then
            ucfj = ucf(j)
            call setcof (eol, j)
            call ypot (j)
            call xpot (j)
            call dacon
            ypj(:n) = yp(:n)
            xpj(:n) = xp(:n)
            xqj(:n) = xq(:n)
            jlast = j
         endif
!
         if (m /= mlast) then
            ucfm = ucf(m)
            call setcof (eol, m)
            call ypot (m)
            call xpot (m)
            call dacon
            ypm(:n) = yp(:n)
            xpm(:n) = xp(:n)
            xqm(:n) = xq(:n)
            mlast = m
         endif
!
         mtp = max(mf(j),mf(m))
!
         if (lfix(m)) then
            ta(1) = 0.d0
            do i = 2, mtp
               ta(i) = rpor(i)*((pf(i,m)*xqj(i)-qf(i,m)*xpj(i))*c+       &
                                (pf(i,m)*pf(i,j)+qf(i,m)*qf(i,j))*ypj(i)) 
            end do

            call quad (result)
! note      ..rinti(m,j,1) is symmetric in (j,m)
            rijm = rinti(m,j,1)
            ecv(li) = (result - rijm)*ucfj

! start dbg
!           WRITE (81,*)'1, RESULT, RIJM, UCFJ, ECV, TA' ! dbg
!           WRITE (81,*)RESULT, RIJM, UCFJ, ECV ! dbg
!           DO i = 1, MTP ! dbg
!              WRITE (81,*) i, TA(i), r(i), rp(i) ! dbg
!           ENDDO ! dbg
! end dbg

         else if (lfix(j)) then
            ta(1) = 0.d0
            do i = 2, mtp
               ta(i) = rpor(i)*((pf(i,j)*xqm(i)-qf(i,j)*xpm(i))*c+      &
                                (pf(i,j)*pf(i,m)+qf(i,j)*qf(i,m))*ypm(i)) 
            end do

!start dbg
!           DO i = 1, MTP
!              WRITE (81,*) i, TA(i)
!               WRITE (83,*) i, r(i), rp(i), rpor(i)
!               WRITE (84,*) i, pf(i,j), qf(i,j)
!               write(85,*) i,ypm(i)
!              write(86,*)i,xpm(i),xqm(i)
!           ENDDO
! end dbg
            call quad (result)

            rijm = rinti(j,m,1)                  !/ nprocs
            ecv(li) = (result - rijm)*ucfm
!start dbg
!           WRITE (81,*)'2, RESULT, RIJM, UCFM, ECV, TA'
!           WRITE (81,*)RESULT, RIJM, UCFJ, ECV, r(i), rp(i)
!end dbg

         else
            
            qdif = abs((ucfj - ucfm)/max(ucfj,ucfm))
            
            if (qdif > p001) then
               obqdif = 1.d0/ucfj - 1.d0/ucfm
               ta(1) = 0.d0
               do i = 2, mtp
                  ta(i) = rpor(i)*((pf(i,m)*xqj(i)-qf(i,m)*xpj(i)          &
                                   -pf(i,j)*xqm(i)+qf(i,j)*xpm(i))*c       &
                          +(ypj(i)-ypm(i))*(pf(i,m)*pf(i,j)+qf(i,m)*qf(i,j))) 
               end do

               call quad (result)
               ecv(li) = result/obqdif

!start dbg
!           WRITE (81,*)'3, RESULT, OBQDIF, ECV, TA'
!           WRITE (81,*)RESULT, OBQDIF, ECV
!           DO i = 1, MTP
!              WRITE (81,*) i, TA(i), r(i), rp(i)
!           ENDDO
!end dbg

            else
 
               obqsum = 1.d0/ucfj + 1.d0/ucfm

               ta(1) = 0.d0
               do i = 2, mtp
                  ta(i) = rpor(i)*((pf(i,m)*xqj(i)-qf(i,m)*xpj(i)         &
                                   +pf(i,j)*xqm(i)-qf(i,j)*xpm(i))*c      &
                         +(ypj(i)+ypm(i))*(pf(i,m)*pf(i,j)+qf(i,m)*qf(i,j))) 
               end do

               call quad (result)
               rijm = rinti(m,j,1)               
               ecv(li) = (result - 2.0*rijm)/obqsum 


!start dbg
!           WRITE (81,*)'4, RESULT, RIUJM, OBQSUM, ECV, TA'
!           WRITE (81,*)RESULT, RIUJM, OBQSUM, ECV
!           DO i = 1, MTP
!              WRITE (81,*) i, TA(i), r(i), rp(i)
!           ENDDO
!end dbg

            endif
         endif

         write (*, 307) np(j), nh(j), np(m), nh(m), ecv(li)

      end do

  302 format(/,'Lagrange multipliers are not required')
  304 format(/,'Include Lagrange multipliers between:'/)
  305 format(13x,2(2x,1i2,1a2))
  306 format(/,'Lagrange multipliers:'/)
  307 format(13x,2(2x,1i2,1a2),2x,1p,d16.9)

      return
      end subroutine setlag
