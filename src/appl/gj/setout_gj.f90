!***********************************************************************
!                                                                      *
      subroutine setout_gj(name, nci, gjc_diag, dgjc_diag)
!                                                                      *
!   Writes diagonal Landé-factor data to <name>.gj or <name>.cgj.     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use def_c,           only: cvac
      use jlabl_c,         only: labj=>jlbr, labp=>jlbp
      use nsmdat_c,        only: sqn, dmomnm, qmomb,                  &
                                 hfsi=>sqn, hfsd=>dmomnm, hfsq=>qmomb
      use prnt_c,          only: nvec, ivec
      use syma_c,          only: iatjpo, iaspar
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=24), intent(in) :: name
      integer, intent(in) :: nci
      real(double), intent(in) :: gjc_diag(nvec), dgjc_diag(nvec)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ierr, jj, k
      real(double) :: fj, gja1, gj, dgj
      character :: filnam*256, form*11, status*3
!-----------------------------------------------
!
      k = index(name,' ')
      if (nci == 0) then
         filnam = name(1:k-1)//'.cgj'
      else
         filnam = name(1:k-1)//'.gj'
      endif
!
      form = 'FORMATTED'
      status = 'NEW'
      call openfl(29, filnam, form, status, ierr)
      if (ierr /= 0) then
         write (6, *) 'Error when opening', filnam
         stop
      endif
!
      write (29, 302) hfsi
      write (29, 303) hfsd
      write (29, 304) hfsq
      write (29, 402)
!
      do i = 1, nvec
         jj = iatjpo(i)
         if (jj <= 1) cycle
         fj = 0.5d00*dble(jj - 1)
         gja1 = sqrt(1.0d00/(fj*(fj + 1.0d00)))
         gj = cvac*gja1*gjc_diag(i)
         dgj = 0.001160d0*gja1*dgjc_diag(i)
!
         write (29, 403) ivec(i), labj(jj), labp((iaspar(i)+3)/2), gj, dgj, &
            gj + dgj
      end do
!
      close(29)
      return
!
  302 format('Nuclear spin                        ',1p,d22.15,' au')
  303 format('Nuclear magnetic dipole moment      ',1p,d22.15,' n.m.')
  304 format('Nuclear electric quadrupole moment  ',1p,d22.15,' barns')
  402 format(/,/,' Interaction constants:'/,/,' Level1  J Parity ',10x,'g_J', &
         14x,'delta g_J',11x,'total g_J'/)
  403 format(1x,1i3,5x,2a4,1p,3d20.10)
      end subroutine setout_gj
