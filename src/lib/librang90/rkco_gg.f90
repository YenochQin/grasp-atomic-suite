module rkco_gg_i
    implicit none (type, external)
    private
    public :: rkco_gg
contains
!***********************************************************************
!                                                                      *
      subroutine rkco_gg (ja,jb,coefficient_consumer,incor,icolbrei)
!                                                                      *
!   Configurations JA, JB. Analyse the tables of quantum numbers set   *
!   in the COMMON  blocks M0 , M1, M2, M3  to determine all possible   *
!   sets of interacting  orbitals which give a non-vanishing Coulomb   *
!   matrix element,  and  initiates the calculation of coefficients.   *
!   The following conventions are in force: (1) labels 1, 2 refer to   *
!   left, right sides of matrix element respectively;   (2) pointers   *
!   JA1, JB1, JA2, JB2 point to the JLIST array of active  orbitals;   *
!   IA1, IB1, IA2, IB2 point to the complete list of orbitals.         *
!                                                                      *
!   Call(s) to: [LIB92]: COR, CORD, ISPAR, ITJPO, SETQNA, VIJOUT.      *
!                                                                      *
!   Rewrite by G. Gaigalas                                             *
!   Transform to fortran 90/95 by G. Gaigalas           December 2012  *
!   The last modification made by G. Gaigalas           October  2020  *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use m_c, only: double, jjc1, jjc2, jjq1, jjq2, jlist, klist, max_orbitals, ncore, npeel, nq1, nq2
      use orb_c, only: double, max_orbitals, orbital_count, orbital_kappa
      use dumx_c, only: double, jc1s, jc2s, jlis, max_orbitals
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE cord_I
      use itjpo_i, only: itjpo
      use ispar_i, only: ispar
      use ichkq2_i, only: ichkq2
      use setqna_i, only: setqna
      use el1_i, only: el1
      use el2_i, only: el2
      use el3_i, only: el3
      use el4_i, only: el4
      use el5_i, only: el5
      use cord_i, only: cord_callback => cord
      implicit none (type, external)
      procedure(cord_callback) :: coefficient_consumer
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      EXTERNAL CORD
      integer, intent(in) :: ja,jb,incor,icolbrei
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i,ib1,ib2,idq,idqg,ii,i1,ij,ijw,im,ipca,it1,it2, &
                 ja1,jb1,ja2,jb2,j,jw,jt1,jt2,jt3,ndq, &
                 klast,kw,kwa,k1,kw1,kw2,npeelm
!-----------------------------------------------
!
!   The Hamiltonian is an even scalar operator
      if (itjpo (ja) /= itjpo (jb)) return
      if (ispar (ja) /= ispar (jb)) return
      if (ichkq2(ja,jb) == 0) return
      call setqna (ja,jb)
!
!   1.0 Analyse peel shell interactions
!
!   1.1 Analyse electron distribution in peel. (The full procedure is
!       needed only if the number of peel orbitals NPEEL .GE. 2)
      if (orbital_count < 1) then
         print *, "RKCO_GG: No subshells."
         stop
      end if
      if (npeel == 0) goto 48
      if (npeel == 1) goto 43
!
!   Find differences in occupations, NDQ, for each peel orbital in
!   turn and use to set up labels of active orbitals maintaining the
!   convention JA1 .LE. JB1, JA2 .LE. JB2.
      idq = 0
      ja1 = 0
      jb1 = 0
      ja2 = 0
      jb2 = 0
      idqg = 0
      do jw = 1,npeel
         j = jlist(jw)
         ndq = nq1(j) - nq2(j)
         if (iabs (ndq) > 2) return
         if (ndq < 0) then
           if (ndq+1 > 0) then
             cycle
           else if (ndq+1 == 0) then
             if (ja2 > 0) then
               jb2 = jw
             else
               ja2 = jw
             end if
             idq = idq+1
             idqg=1+idqg
           else if (ndq+1 < 0) then
             ja2 = jw
             idq = idq+2
             idqg=20+idqg
           end if
         else
           if (ndq-1 > 0) then
             ja1 = jw
             idq = idq+2
             idqg=20+idqg
             cycle
           else if (ndq-1 == 0) then
             if (ja1 > 0) then
               jb1 = jw
             else
               ja1 = jw
             end if
             idq = idq+1
             idqg=1+idqg
             cycle
           else
             cycle
           end if
         end if
      end do
!
!   1.2 Calculate coefficients for all possible sets of active shells.
!
!   There are 4 cases, depending on the value of IDQ, the sum of the
!   absolute differences NDQ:
!
!   1.2.1 IDQ .GT. 4: matrix element null
      if (idq > 4) return
      if (idq == 4) then
!
!   1.2.2 IDQ .EQ. 4: matrix element uniquely defined
        if (jb1 == 0) then
          jb1 = ja1
        end if
        if (jb2 == 0) then
          jb2 = ja2
        end if
        continue
        if(idqg/=40) then
          if(idqg/=22) then
!
!         TARP KONFIGURACIJU
!                        KAI      N'1=N1+-1
!                                 N'2=N2+-1
!                        KAI      N'3=N3+-1
!                                 N'4=N4+-1
            call el5(ja,jb,ja1,jb1,ja2,jb2,icolbrei)
          else
            call el4(ja,jb,ja1,jb1,ja2,jb2,icolbrei)
          end if
        else
!
!       TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
          call el2(ja,jb,ja1,ja2,icolbrei)
        end if
        return
      end if
!
      if (idq /= 2) then
        if (idq /= 0) then
!
!         3.0 Diagnostic print - NW .LT. 1
          write (*,300)
          stop
        end if
      else
        klast = 1
        goto 16
      end if
!
      if (ja == jb) goto 43
      klast = npeel
!
!   1.2.3 IDQ .EQ. 2: One orbital fixed each side include all
!                     possible spectators.
!
!   Also IDQ .EQ. 0 for a matrix element off-diagonal in coupling
!   only. Must sum over all pairs of orbitals excluding core-core
!   terms
   16 do kwa = 1,klast
         if (idq /= 2) then
           ja1 = kwa
           ja2 = kwa
         end if
         jt1 = ja1
         jt2 = ja2
         it1 = jlist(ja1)
         it2 = jlist(ja2)
         do kw = kwa,npeel
            k1 = jlist(kw)
            if (nq1(k1)*nq2(k1) == 0) cycle
            jb1 = kw
            jb2 = kw
            ja1 = jt1
            ja2 = jt2
!
!   Interchange JA1 and JB1 and/or JA2 and JB2 if necessary
            if (ja1-jb1 > 0) then
              jt3 = jb1
              jb1 = ja1
              ja1 = jt3
            else if (ja1-jb1 == 0) then
              ib1 = jlist(jb1)
              if (nq1(ib1) <= 1) cycle
            end if
            if (ja2-jb2 > 0) then
              jt3 = jb2
              jb2 = ja2
              ja2 = jt3
            else if (ja2-jb2 == 0) then
              ib2 = jlist(jb2)
              if (nq2(ib2) <= 1) cycle
            end if
            if(idq/=0) then
!
!     TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
!                                 N'3=N3
              call el3(ja,jb,ja1,jb1,ja2,jb2,icolbrei)
            else
!
!     TARP TU PACIU KONFIGURACIJU
              call el1(ja,jb,ja1,jb1,1,icolbrei)
            end if
         end do
         if ((idq == 0) .and. (ncore == 0)) cycle
         if ((ncore == 0) .or. (orbital_kappa(it1) /= orbital_kappa(it2))) return
!
!   This section calculates the terms arising from active electrons
!   which are in closed shells
         npeelm = npeel-1
         do i = 1,npeel
            jlis(i) = jlist(i)
         end do
         do i = 1,npeelm
            jc1s(i) = jjc1(i)
            jc2s(i) = jjc2(i)
         end do
         do kw = 1,ncore
            ijw = klist(kw)
            do i = 1,npeel
               ij = jlist(i)
               if (ijw < ij) goto 29
            end do
            i = npeel+1
            goto 31
   29       im = npeel-i+1
            do ii = 1,im
               jlist(npeel+2-ii) = jlist(npeel+1-ii)
               if (npeel == ii) goto 31
               jjc1(npeel+1-ii) = jjc1(npeel-ii)
               jjc2(npeel+1-ii) = jjc2(npeel-ii)
            end do
   31       continue
            if (i < 3) then
              i1 = jlist(1)
              jjc1(1) = jjq1(3,i1)
              jjc2(1) = jjq2(3,i1)
            else
              jjc1(i-1) = jjc1(i-2)
              jjc2(i-1) = jjc2(i-2)
            end if
            jlist(i) = ijw
            ja1 = jt1
            if (jt1 >= i) ja1 = ja1+1
            jb1 = i
            ja2 = jt2
            if (jt2 >= i) ja2 = ja2+1
            jb2 = i
            if (ja1-jb1 > 0) then
              jt3 = jb1
              jb1 = ja1
              ja1 = jt3
            end if
            if (ja2-jb2 > 0) then
              jt3 = jb2
              jb2 = ja2
              ja2 = jt3
            end if
            npeel = npeel+1
            if(idq/=0) then
              if(idqg/=40) then
                if(idqg/=2) then
                   write(99,995)
                   return
                else
!
!     TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
!                                 N'3=N3
                  call el3(ja,jb,ja1,jb1,ja2,jb2,icolbrei)
                end if
              else
                write(99,994)
                return
              end if
            else
!
!     TARP TU PACIU KONFIGURACIJU
               call el1(ja,jb,ja1,jb1,1,icolbrei)
            end if
            npeel = npeel-1
            npeelm = npeel-1
            do i = 1,npeel
               jlist(i) = jlis(i)
            end do
            do i = 1,npeelm
               jjc1(i)  = jc1s(i)
               jjc2(i)  = jc2s(i)
            end do
         end do
      end do
      return
!
!   1.2.4 IDQ .EQ. 0 - diagonal case. Include all pairs with
!         JA1 = JA2, JB1 = JB2.
   43 do kw1 = 1,npeel
         k1 = jlist(kw1)
         jb1 = kw1
         jb2 = kw1
         do kw2 = 1,kw1
            ja1 = kw2
            if (ja1 == jb1) then
              if (nq1(k1) <= 1) cycle
            end if
            ja2 = ja1
            if(ja /= jb) then
              if(npeel == 1 .and. nq1(k1) == nq2(k1)) then
!
!               TARP TU PACIU BUSENU
                call el1(ja,jb,ja1,jb1,1,icolbrei)
              else if(idqg/=2) then
                write(99,996)
                write(99,*)"JA,JB,JA1,JB1",ja,jb,ja1,jb1
                write(99,*)"IDQG IDQ",idqg,idq
                return
              else
!
!                TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
!                                 N'3=N3
                 call el3(ja,jb,ja1,jb1,ja2,jb2,icolbrei)
              end if
            else
!
!             TARP TU PACIU BUSENU
              call el1(ja,jb,ja1,jb1,0,icolbrei)
            end if
         end do
      end do
   48 if (incor < 1) return
      if (ncore == 0) return
!
!   2.0 The diagonal case. deal with contributions from core orbitals
!       if INCOR .EQ. 1.
      if(ja/=jb) return
      do kw1 = 1,ncore
         jb1 = kw1
         jb2 = kw1
!
!   2.1 Calculate contribution from core/core terms
         ipca = 2
         do kw2 = 1,kw1
            ja1 = kw2
            ja2 = kw2
            call coefficient_consumer (ja,jb,ja1,ipca,jb1)
         end do
!
!   2.2 Calculate contribution from peel/core terms
         if (npeel == 0) cycle
         ipca = 1
         do kw2 = 1,npeel
            ja1 = kw2
            ja2 = kw2
            call coefficient_consumer (ja,jb,ja1,ipca,jb1)
         end do
      end do
      return
  300 format ("RKCO_GG: Error.")
  994 format("   rie zymes 38?? atv N=N-N  !!!!!!")
  995 format("   rie zymes 38 atv N=N-N  !!!!!!")
  996 format("   rie zymes 45 atv N=N-N  !!!!!!")
      end subroutine rkco_gg
end module rkco_gg_i
