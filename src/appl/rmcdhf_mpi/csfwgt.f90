!***********************************************************************
      subroutine csfwgt(lstdio)
!                                                                      *
!   Print  the  weights of the largest five CSF contributors to each   *
!   ASF.                                                               *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, ISPAR, ITJPO.          *
!                                                                      *
!                                          Last updated: 21 Dec 1992   *
!                                          Last updated: 24 Feb 1997   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use def_c
      use eigv_c
      use jlabl_c, labj=>jlbr, labp=>jlbp
      use syma_c,          only: jpgg
      use iounit_c
      use hblock_c
      use pos_c
      use blkidx_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) :: lstdio
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(5) :: iconf
      integer :: iv, jblock, ncf, ncfpat, ncminpat, nevecpat, nevecoff, nelt, &
         icf, ivtjpo, ivspar, i, j, itemp, irem, ip
      real(double) :: temp, w
      real(double), dimension(5) :: coeff
      character :: record*256, cnum*8
!-----------------------------------------------
!

      if (lstdio) then
         write (istdo, 300)
      else
         write (24, 300)
      endif

      do iv = 1, ncmin
!         loop over eigenvectors

         jblock = idxblk(iv)
         ncf = ncfblk(jblock)
         ncfpat = ncfpast(jblock)
         ncminpat = ncminpast(jblock)
         nevecpat = nevecpast(jblock)
         nevecoff = nevecpat + (iv - ncminpat - 1)*ncf
         nelt = min(5,ncf)                       !Find maximum 5 ... within block
         icf = iccmin(iv)
!GGGG
         ivtjpo = iabs(jpgg(jblock))
         if(jpgg(jblock) .ge. 0) then
            ivspar = 2
         else
            ivspar = 1
         end if
!GG         IVTJPO = IATJPO(IV)                     ! j-value related
!GG         IVSPAR = IASPAR(IV)                     ! parity related


         do i = 1, nelt
            coeff(i) = evec(nevecoff + i)
            iconf(i) = i
         end do

!         sort the first nelt in decreasing order
         do i = 1, nelt
            do j = i + 1, nelt
               if (abs(coeff(j)) <= abs(coeff(i))) cycle
               temp = coeff(i)
               coeff(i) = coeff(j)
               coeff(j) = temp
               itemp = iconf(i)
               iconf(i) = iconf(j)
               iconf(j) = itemp
            end do
         end do

         l20: do i = nelt + 1, ncf
            w = evec(nevecoff + i)
            if (w==0.d0 .or. abs(w)<=abs(coeff(nelt))) cycle  l20
!               we have a non-zero value larger than the largest so far
            do j = 1, nelt
               if (abs(w) <= abs(coeff(j))) cycle
               coeff(nelt:1+j:(-1)) = coeff(nelt-1:j:(-1))
               iconf(nelt:1+j:(-1)) = iconf(nelt-1:j:(-1))
               coeff(j) = w
               iconf(j) = i
               cycle  l20
            end do
         end do l20

!GG         IP = (IASPAR(IV) + 3)/2
           ip =  ivspar

         if (lstdio) then
            write (istdo, 320) jblock, icf, labj(ivtjpo), labp(ip), (coeff(i),i&
               =1,nelt)
            write (istdo, 330) (iconf(i),i=1,nelt)
         else

            write (24, 320) jblock, icf, labj(ivtjpo), labp(ip), (coeff(i),i=1,&
               nelt)
            write (24, 330) (iconf(i),i=1,nelt)
         endif
      end do

  300 format(/,'Weights of major contributors to ASF:'/,/,&
         'Block Level J Parity      CSF contributions'/)
  310 format(1x,a14,80a)
  320 format(i3,1x,i5,2x,2a4,5(3x,f8.4))
  330 format(19x,5(3x,i8))

      return
      end subroutine csfwgt
