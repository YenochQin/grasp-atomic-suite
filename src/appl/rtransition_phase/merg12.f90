!***********************************************************************
!                                                                      *
      subroutine merg12(name, ncorer, ncore)
!                                                                      *
!   This subroutines merges the initial and final state lists          *
!   Observ that there may doublets in this list if the initial         *
!   and final states have the same parity                              *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,   only: max_orbitals
      use def_c
      use orb_c
      use blk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i
      use iqr_i
      use ispar_i
      use isparr_i
      use itjpo_i
      use itjpor_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ncorer
      integer , intent(in) :: ncore
      character , intent(in) :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j1
      integer :: newnp(0:max_orbitals + 1)
      integer :: newnak(max_orbitals), icon(max_orbitals), icon1(max_orbitals)
      integer :: i, j, icomp, new, ipi1, ipi2, k, new1, new2, m, nline
      character :: line*500
      character :: newnh(max_orbitals)*2
      character :: new3*2
!-----------------------------------------------
!
!
! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
! NCFF(I): the end position of the Ith block for the final states in the globle CSF list

!
      open(unit=21, file='SLASK', form='FORMATTED', status='UNKNOWN', position=&
         'asis')
!
!   The same number of electrons must appear in both lists
!
      if (nelecr /= electron_count) then
         write (6, *) 'The number of electrons is not equal in the'
         write (6, *) ' first and second GRASP92 Configuration'
         write (6, *) ' symmetry list files.'
         stop
      endif
!
!   THe core orbitals must be the same
!
      if (ncore /= ncorer) then
         write (6, *) 'The number of core orbitals must be the same'
         stop
      endif

      do i = 1, ncore
         if (orbital_principal_n(i)==npr(i) .and. orbital_kappa(i)==nakr(i)) cycle
         write (6, *) 'The core orbitals must be the same'
         stop
      end do

      icon(:orbital_count) = 0
!
!   For each orbital in the initial list check if there is a corresponding
!   orbital in the final list. If so give the number of the corresponding
!   orbital
!
      do i = 1, orbital_count
         do j = 1, nwr
            if (orbital_principal_n(i)/=npr(j) .or. orbital_kappa(i)/=nakr(j)) cycle
            icon(i) = j
         end do
      end do
!
!   Check if the ordering of the initial and final state orbitals
!   is consistent. The condition for this is that ICON is increasing
!
      j = 0
      do i = 1, orbital_count
         if (icon(i) == 0) cycle
         j = j + 1
         icon1(j) = icon(i)
      end do

      icomp = icon1(1)
      do i = 2, j
         if (icon1(i) < icomp) then
            write (*, *) ' In merg12: ordering of the initial and final'
            write (*, *) ' state orbitals is inconsistent. STOP'
            stop
         else
            icomp = icon1(i)
         endif
      end do
!
!   Determine a common orbital set for the initial and final state.
!   The common set must be such that the order of both the initial and
!   final state sets are preserved.
!
      do i = 1, nwr
         newnp(i) = npr(i)
         newnak(i) = nakr(i)
         newnh(i) = nhr(i)
         write (*, *) newnp(i), newnh(i)
      end do
!
!   Add the initial state orbitals at the end
!
      new = nwr
      do i = 1, orbital_count
         if (icon(i) /= 0) cycle
         new = new + 1
         newnp(new) = orbital_principal_n(i)
         newnak(new) = orbital_kappa(i)
         newnh(new) = orbital_symmetry_labels(i)
         write (*, *) newnp(new), newnh(new), new
      end do
!
!   Now sort in the orbitals at the end in the right position
!
      l193: do i = nwr + 1, new
!
!   Position in initial state list
!
         do j = 1, orbital_count
            if (newnp(i)/=orbital_principal_n(j) .or. newnak(i)/=orbital_kappa(j)) cycle
            ipi1 = j
         end do
         write (*, *) 'i,ipi1', i, ipi1

         do j = 1, nwr
            ipi2 = 0
            do k = 1, orbital_count
               if (newnp(j)/=orbital_principal_n(k) .or. newnak(j)/=orbital_kappa(k)) cycle
               ipi2 = k
            end do
            write (*, *) 'j,ipi2', j, ipi2
            if (ipi2 == 0) cycle
            if (ipi1 >= ipi2) cycle
            new1 = newnp(i)
            new2 = newnak(i)
            new3 = newnh(i)
            newnp(i:1+j:(-1)) = newnp(i-1:j:(-1))
            newnak(i:1+j:(-1)) = newnak(i-1:j:(-1))
            newnh(i:1+j:(-1)) = newnh(i-1:j:(-1))
            newnp(j) = new1
            newnak(j) = new2
            newnh(j) = new3
            do m = 1, new
               write (*, *) newnp(m), newnh(m)
            end do
            cycle  l193
         end do
      end do l193

      orbital_count = new
      orbital_principal_n(:orbital_count) = newnp(1:orbital_count)
      orbital_kappa(:orbital_count) = newnak(:orbital_count)
      orbital_symmetry_labels(:orbital_count) = newnh(:orbital_count)
!
!   Determine NKL and NKJ
!
      do i = 1, orbital_count
         orbital_two_j(i) = 2*abs(orbital_kappa(i)) - 1
         if (orbital_kappa(i) > 0) then
            orbital_l(i) = (orbital_two_j(i)+1)/2
         else
            orbital_l(i) = (orbital_two_j(i)-1)/2
         endif
      end do
!
!   Determine the common core subshells; write out the list;
!   determine the pell subshells; write out the list; these
!   data form the first part of the header of the .csl file;
!   one additional line forms the remainder of the header of
!   the  .csl file
!
      write (21, '(A)') 'Core subshells:'
      write (21, 301) (orbital_principal_n(i),orbital_symmetry_labels(i),i=1,ncore)
      write (21, '(A)') 'Peel subshells:'
      write (21, 301) (orbital_principal_n(i),orbital_symmetry_labels(i),i=ncore + 1,orbital_count)
      write (21, '(A)') 'CSF(s):'
!
!   Now write out all CSFs in the initial and final state list
!
      j = index(name(1),' ')
      open(unit=23, file=name(1)(1:j-1)//'.c', form='FORMATTED', status='OLD', &
         position='asis')

      do i = 1, 5
         read (23, '(A)') line
      end do

      nblocki = 0
      nline = 0
    5 continue
      read (23, '(A)', end=98) line
      if (line(1:2) == ' *') then
         nblocki = nblocki + 1
         ncfi(nblocki) = nline/3
      else
         nline = nline + 1
      endif
      k = 500
   10 continue
      if (line(k:k) == ' ') then
         k = k - 1
         if (k > 1) go to 10
      endif
      write (21, '(A)') line(1:k)
      go to 5
   98 continue
      nblocki = nblocki + 1
      ncfi(nblocki) = nline/3
      close(23)

! zou
!     if(NAME(2).EQ.NAME(1)) return
! zou
      j = index(name(2),' ')
      open(unit=23, file=name(2)(1:j-1)//'.c', form='FORMATTED', status='OLD', &
         position='asis')

      do i = 1, 5
         read (23, '(A)') line
      end do

      nblockf = 0
      nline = 0
   15 continue
      read (23, '(A)', end=99) line
      if (line(1:2) == ' *') then
         nblockf = nblockf + 1
         ncff(nblockf) = nline/3
      else
         nline = nline + 1
      endif
      k = 500
   20 continue
      if (line(k:k) == ' ') then
         k = k - 1
         if (k > 1) go to 20
      endif
      write (21, '(A)') line(1:k)
      go to 15
   99 continue
      nblockf = nblockf + 1
      ncff(nblockf) = nline/3
      close(23)

      close(21)
!

      write (6, *) nblocki
      write (6, *) (ncfi(i),i=1,nblocki)
      write (6, *) nblockf
      write (6, *) (ncff(i),i=1,nblockf)
  301 format(120(1x,1i2,1a2))
!
      return
      end subroutine merg12
