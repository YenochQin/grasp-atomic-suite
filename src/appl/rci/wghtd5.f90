!***********************************************************************
!                                                                      *
      subroutine wghtd5(iatjpo, iaspar)
!                                                                      *
!   Print  the  weights of the largest five CSF contributors to each   *
!   ASF.                                                               *
!                                                                      *
!   Call(s) to: ALLOC, DALLOC, ISPAR, ITJPO.                           *
!                                                                      *
!                                          Last updated: 02 Nov 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use eigv_c
      use jlabl_c, labj=>jlbr, labp=>jlbp
      use orb_c, only: csf_count, orbital_count, csf_subshell_occupations
      use prnt_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE ispar_I
!      USE itjpo_I
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: iatjpo, iaspar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer, dimension(5) :: iconf
!      REAL(DOUBLE), DIMENSION(5) :: WGHT(5)
      real(double), dimension(5) :: wt(5)
!      REAL(DOUBLE), DIMENSION(:), pointer :: wt
!      INTEGER, DIMENSION(:), pointer :: next
      real(double) :: awt
      integer :: i, j, jj, ibegin, iend, nelt, iv, icf, ip
!
!   Allocate storage for local arrays
!
!CFF  These arrays are not needed
!      CALL ALLOC (WT, NCF,'WT', 'WGHTD5')
!      CALL ALLOC (NEXT,NCF,'NEXT', 'WGHTD5')
!
      write (24,300)
!
!CFF  Initialized NELT for NCF >5
      nelt = 5
      if (csf_count .lt. 5) nelt = csf_count
      do 7 iv = 1,nvec
!
         icf = ivec(iv)
!        .., initialize for a new eigenvector
         wt = 0.d0;  iconf = 0
         ibegin = (iv-1)*csf_count
         iend   = iv*csf_count
!
         do 4 i = ibegin+1, iend
            awt = abs(csf_mixing_coefficients(i))
            if (awt > wt(5)) then
!              ... the weight needs to be inserted
               do j = 1,5
                  if ( awt > wt(j) ) then
                     if (j < 5) then
!                       ... shift right
                        do jj = 5, j+1,-1
                           wt(jj) = wt(jj-1)
                           iconf(jj) = iconf(jj-1)
                        end do
                     end if
!                    ...  insert the element
                     wt(j) = awt
                     iconf(j) = i
                     exit
                  end if
               end do
            end if
    4    continue
!
!   Print first five elements of list.
!
         ip = (iaspar + 3 ) / 2
!        NELT = MIN (I,5)
         write (24,301) icf,labj(iatjpo),labp(ip),      &
               (csf_mixing_coefficients(iconf(j)),j = 1,nelt)
         write (24,302) (iconf(j)-ibegin,j = 1,nelt)
    7 continue
!
!   Deallocate storage for local arrays
!
!      CALL DALLOC (WT, 'WT', 'WGHTD5')
!      CALL DALLOC (NEXT, 'NEXT', 'WGHTD5')
!
      return
!
  300 format (/'Weights of major contributors to ASF:'                  &
             //'Level J Parity      CSF contributions'/)
  301 format (i3,2x,2a4,5(3x,f8.5))
  302 format (13x,5(i11))
!
      end subroutine wghtd5
