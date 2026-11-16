module oneparticlejj_i
    implicit none (type, external)
    private
    public :: oneparticlejj
contains
!*******************************************************************
!                                                                  *
      subroutine oneparticlejj(ka,iopar,ja,jb,ia1,ia2,vshell)
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
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only:  max_orbitals
      use cons_c
      use m_c
      use orb_c,   only: orbital_count, orbital_kappa
      use dumx_c,  only: jlis, jc1s, jc2s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use itrig_i, only: itrig
      use ichkq1_i, only: ichkq1
      use ichop_i, only: ichop
      use ispar_i, only: ispar
      use itjpo_i, only: itjpo
      use setqna_i, only: setqna
      use oneparticlejj1_i, only: oneparticlejj1
      use oneparticlejj2_i, only: oneparticlejj2
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: ka,iopar,ja,jb
      integer, intent(out) :: ia1,ia2
      real(real64), intent(out) :: vshell(max_orbitals)
!      DIMENSION VSHELL(NNNW)
!      DIMENSION IS(2),KS(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kk,ij,idq,ja1,ja2,jw,ndq,ns,ish,i,ii,i1,im,  &
                 ks1,ks2,nx,npeelm
      integer :: is(2),ks(2)
      real(real64) :: tcoeff
!-----------------------------------------------
      ia1 = 0
      kk = ka+ka+1
      if(itrig(itjpo(ja),itjpo(jb),kk)==0)return
      if((iopar/=0).and.(ispar(ja)*ispar(jb)*iopar/=1))return
      if(ichkq1(ja,jb)==0)return
!
      call setqna (ja,jb)
!
      do ij = 1,orbital_count
         vshell(ij) = zero
      end do
!
!   Analyse peel shell interactions
      idq = 0
      ja1 = 0
      ja2 = 0
      if (npeel /= 0) then
        do jw = 1,npeel
          ij = jlist(jw)
          ndq = nq1(ij)-nq2(ij)
          if (abs (ndq) > 1) return
          if (ndq > 0) then
            ja1 = jw
            idq = idq+1
          else if (ndq < 0) then
            ja2 = jw
            idq = idq+1
          end if
        end do
        if (idq > 2) return
!
!   Evaluate the array VSHELL
!
!   Then there are two possibilities IDQ = 0 or IDQ = 2
!   if IDQ = 0, then loop over all shells by index ISH
!   if IDQ = 2, then one orbital fixed on each side
        ns = npeel
      end if
      if (idq == 2) goto 19
!
!   Loop over shells when IDQ = 0
      ish = 0
      if (npeel == 0) goto 9
      do i = 1,npeel
         jlis(i) = jlist(i)
      end do
      if (npeel == 1) goto 9
      npeelm = npeel-1
      do i = 1,npeelm
         jc1s(i) = jjc1(i)
         jc2s(i) = jjc2(i)
      end do
!
!   If ISH .GT. NW, then loop is over and return
    9 ish = ish+1
      if (ish > orbital_count) return
      if (ichop (ish,ja) == -1) goto 9
      if (ichop (ish,ja) == 0) goto 16
!
!   Case one --- the ISH-th shell is in the core or in the peel and
!   closed for both sides
      i = 1
      if (npeel==0) goto 15
      do i = 1,npeel
        ij = jlist(i)
        if (ish < ij) goto 11
      end do
      i = npeel+1
      goto 13
   11 im = npeel-i+1
      do ii = 1,im
         jlist(npeel+2-ii) = jlist(npeel+1-ii)
         if (npeel==ii) goto 13
         jjc1(npeel+1-ii) = jjc1(npeel-ii)
         jjc2(npeel+1-ii) = jjc2(npeel-ii)
      end do
   13 continue
      if (i < 3) goto 14
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
      do jw = 1,npeel
        nx = ish-jlist(jw)
        if (nx==0) goto 18
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
      ks1 = 2*abs (orbital_kappa(ia1))
      ks2 = 2*abs (orbital_kappa(ia2))
!
!   Check triangular condition for the active shells
      if (itrig (ks1,ks2,kk)==1) goto 99
      if (idq == 2) return
      goto 100
!
!   Set tables of quantum numbers of non-interacting spectator shells
   99 continue
      if (idq == 0) then
!GG
!GG   Cia orginalioje programoje yra klaida
!GG
!GG        IF(JA .EQ. JB) THEN
          call oneparticlejj1(ns,ka,ja,jb,ja1,ja2,tcoeff)
!GG          CALL ONESCALAR1(NS,JA,JB,JA1,JA2,TCOEFF)
!GG        ELSE
!GG          TCOEFF = 0.0D 00
!GG        END IF
      else if (idq == 2) then
!
!   IDQ = 2 Case
!
!       Permutation factor for IDQ = 2
        call oneparticlejj2(ns,ka,ja1,ja2,tcoeff)
!GG        CALL ONESCALAR2(JA,JB,JA1,JA2,TCOEFF)
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
      if (npeel == 0) goto 9
      do i = 1,npeel
         jlist(i) = jlis(i)
      end do
      if (npeel == 1) goto 9
      npeelm = npeel-1
      do  i = 1,npeelm
         jjc1(i)  = jc1s(i)
         jjc2(i)  = jc2s(i)
      end do
      goto 9
      return
      end subroutine oneparticlejj
end module oneparticlejj_i
