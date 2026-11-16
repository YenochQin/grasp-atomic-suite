!*******************************************************************
!                                                                  *
      subroutine onescalar(ja,jb,ia1,ia2,vshell)
!                                                                  *
!   The main program for evaluating the reduced matrix elements of *
!   a one particle operator for configurations in jj-coupling.     *
!                                                                  *
!   Call(s) to: [LIB92]: CFP, FIXJ, GENSUM, ICHOP, IROW1, ISPAR,   *
!                        ITJPO, ITRIG, SETQNA, VIJOUT.             *
!               [NJGRAF]: NJGRAF.                                  *
!                                                                  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use parameter_def,   only:  nnnw
      use cons_c
      use m_c
      use orb_c,   only: nw, nak
      use dumx_c,  only: jlis, jc1s, jc2s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use itrig_i
      use ichkq1_i
      use ichop_i
      use ispar_i
      use itjpo_i
      use setqna_i
      use onescalar1_i
      use onescalar2_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ja,jb
      integer, intent(out)      :: ia1,ia2
      real(double), intent(out) :: vshell(nnnw)
!      DIMENSION VSHELL(NNNW)
!      DIMENSION IS(2),KS(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kk,iopar,ij,idq,ja1,ja2,ndq,ns,ish,i,ii,i1,im, &
                 jw,ks1,ks2,nx,npeelm
      integer, dimension(2) :: is,ks
      real(double) ::  tcoeff
!-----------------------------------------------
      ia1 = 0
      kk = 1
      iopar = 1
      if(itrig(itjpo (ja),itjpo (jb),kk).eq.0)return
      if((iopar.ne.0).and.(ispar(ja)*ispar(jb)*iopar.ne.1))return
      if(ichkq1(ja,jb).eq.0)return
!
      call setqna (ja,jb)
!
      do ij = 1,nw
         vshell(ij) = zero
      end do
!
!   Analyse peel shell interactions
      idq = 0
      ja1 = 0
      ja2 = 0
      if (npeel .ne. 0) then
        do jw = 1,npeel
          ij = jlist(jw)
          ndq = nq1(ij)-nq2(ij)
          if (iabs (ndq) .gt. 1) return
          if (ndq .gt. 0) then
            ja1 = jw
            idq = idq+1
          elseif (ndq .lt. 0) then
            ja2 = jw
            idq = idq+1
          endif
        end do
!
        if (idq .gt. 2) return
!
!   Evaluate the array VSHELL
!
!   Then there are two possibilities IDQ = 0 or IDQ = 2
!   if IDQ = 0, then loop over all shells by index ISH
!   if IDQ = 2, then one orbital fixed on each side
        ns = npeel
      endif
!
      if (idq .eq. 2) goto 19
!
!   Loop over shells when IDQ = 0
      ish = 0
      if (npeel .eq. 0) goto 9
      do i = 1,npeel
         jlis(i) = jlist(i)
      end do
      if (npeel .eq. 1) goto 9
      npeelm = npeel-1
      do i = 1,npeelm
         jc1s(i) = jjc1(i)
         jc2s(i) = jjc2(i)
      end do
!
!   If ISH .GT. NW, then loop is over and return
    9 ish = ish+1
      if (ish .gt. nw) return
      if (ichop (ish,ja) .eq. -1) goto 9
      if (ichop (ish,ja) .eq. 0) goto 16
!
!   Case one --- the ISH-th shell is in the core or in the peel and
!   closed for both sides
      i = 1
      if (npeel.eq.0) goto 15
      do i = 1,npeel
        ij = jlist(i)
        if (ish .lt. ij) goto 11
      end do
      i = npeel+1
      goto 13
   11 im = npeel-i+1
      do ii = 1,im
         jlist(npeel+2-ii) = jlist(npeel+1-ii)
         if (npeel.eq.ii) goto 13
         jjc1(npeel+1-ii) = jjc1(npeel-ii)
         jjc2(npeel+1-ii) = jjc2(npeel-ii)
      end do
   13 continue
      if (i .lt. 3) goto 14
      jjc1(i-1) = jjc1(i-2)
      jjc2(i-1) = jjc2(i-2)
      goto 15
   14 i1 = jlist(1)
      jjc1(1) = jjq1(3,i1)
      jjc2(1) = jjq2(3,i1)
   15 jlist(i) = ish
      ja1 = i
      ja2 = i
      ns = npeel+1
      goto 19
!
!   Case two --- the ISH-th shell is in the peel and open for either
!   side
   16 ns = npeel
      do  jw = 1,npeel
        nx = ish-jlist(jw)
        if (nx.eq.0) goto 18
      end do
   18 ja1 = jw
      ja2 = jw
!
!   Main computation
!
!     JA1, JA2 are the indices of interacting shells in JLIST
!     IA1, IA2 are the indices of interacting shells in NW
   19 ia1 = jlist(ja1)
      ia2 = jlist(ja2)
      ks1 = 2*iabs (nak(ia1))
      ks2 = 2*iabs (nak(ia2))
!
!   Check triangular condition for the active shells
      if (itrig (ks1,ks2,kk).eq.1) goto 99
      if (idq .eq. 2) return
      goto 100
!
!   Set tables of quantum numbers of non-interacting spectator shells
   99 continue

      if (idq .eq. 0) then
!GG
!GG   Cia orginalioje programoje yra klaida
!GG
        if(ja .eq. jb) then
          call onescalar1(ns,ja,jb,ja1,ja2,tcoeff)
        else
          tcoeff = zero
        end if
      else if (idq .eq. 2) then
!
!   IDQ = 2 Case
!
!       Permutation factor for IDQ = 2
        call onescalar2(ja,jb,ja1,ja2,tcoeff)
        vshell(1) = tcoeff
        return
      end if
!
!   End of loop over parent states
!
!
!   IDQ = 0 CASE
      vshell(ish) = tcoeff
!
!   Loop over all shells when IDQ = 0
100   continue
      if (npeel .eq. 0) goto 9
      do i = 1,npeel
         jlist(i) = jlis(i)
      end do
      if (npeel .eq. 1) goto 9
      npeelm = npeel-1
      do  i = 1,npeelm
         jjc1(i)  = jc1s(i)
         jjc2(i)  = jc2s(i)
      end do
      goto 9
      return
      end subroutine onescalar
