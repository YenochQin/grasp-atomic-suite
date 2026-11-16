!***********************************************************************
!                                                                      *
      subroutine lodcsh(nfile, ncore)
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Call(s) to: [LIB92]:  PRSRCN, PRSRSL       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by C. F. Fischer to read only the header information
!
! Input:
!   nfile
!
! Output:
!   ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,   only: max_orbitals
      use def_c
      use orb_c
      use terms_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use prsrsl_i
      use prsrcn_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nfile
      integer, intent(out) :: ncore
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iocc(max_orbitals)
      integer :: iqadum, npeel, i, j, npj, nakj, ios, ierr
      character :: str*256
!-----------------------------------------------
!
!   Entry message
!
      write (6, *) 'Loading CSF file ... Header only'
!
!   Get the list of subshells
!
      orbital_count = 0
!
!   Read the list of core subshells; set up the arrays NP, NAK,
!   NKL, NKJ, NH for these subshells
!
      call prsrsl (nfile, 1)
      ncore = orbital_count
!
!   Skip the peel subshell identification header; read the list of
!   peel subshells; set up the arrays NP, NAK, NKL, NKJ, NH for
!   these subshells
!
      read (nfile, *)
      call prsrsl (nfile, 2)
      npeel = orbital_count - ncore
!
!   Ensure that the sets of core and peel subshell are disjoint
!
      do j = ncore + 1, orbital_count
         npj = orbital_principal_n(j)
         nakj = orbital_kappa(j)
         do i = 1, ncore
            if (orbital_principal_n(i)/=npj .or. orbital_kappa(i)/=nakj) cycle
            write (stderr_unit, *) 'lodcsh: The lists of core and', &
               ' peel subshells must form disjoint sets.'
            stop
         end do
      end do

      write (6, *) 'There are/is ', orbital_count, ' relativistic subshells;'
!
!   Skip the header for the list of CSFs
!
      read (nfile, *)
!
!   To determine the number of electrons. This was done very much later
!   in the non-block mcp program( near the end of lodcsh). In the block
!   version, that will be used as a check to the value obtained below.
!
      read (nfile, '(A)', iostat=ios) str
      call prsrcn (str, ncore, iocc, ierr)
      backspace (nfile)                          ! return to the first CSF item
!             Number of electrons in the peel shells
      electron_count = sum(iocc(ncore+1:orbital_count))
!             Add the number of electrons in the core shells
      electron_count = electron_count + sum(orbital_two_j(:ncore)+1)
      return
      end subroutine lodcsh
