!***********************************************************************
!                                                                      *
      subroutine mctin(iopar, jkp, name)
!                                                                      *
!   This routine loads  coefficients with parity and  rank specified   *
!   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
!   (+/- 1) and is determined from the sign of  KP(JKP).               *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:26:50   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix, max_orbitals
      use debug_c
      use decide_c
      use def_c
      use foparm_c
      use offd_c
      use orb_c
      use osc_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcnsa_i
      use alcnta_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: iopar
      integer, intent(in) :: jkp
      character  :: name(2)*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nfile = 93
      integer, parameter :: nfile1 = 237
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nfile2,m,k,iblki,iblkf,i,labl,ncsf,j,nprocsd,myidd
      real(real64) :: tshell(max_orbitals)
      logical :: avail
!-----------------------------------------------
!
!   Read the data back as required by OSCL conventions
!
      nfile2 = nfile1 + jkp

      m = 0
      k = 0
!
      read (nfile2) iblki,iblkf,orbital_count,nkp,nprocsd,myidd
      read (nfile2) nint
!
      do i = 1, nint
!
         read (nfile2) labl, ncsf
!
         m = m + 1
!bieron
         if (m >= nsdim) call alcnsa (jja, jjb, hb1, hb2, hc1, &
            hc2, hm1, hm2, lab, nptr, nsdim, 2)
         lab(m) = labl
!
!   Read configuration pairs and coefficients for this integral
!
    4    continue
         if (ncsf + k > ntdim) then
            call alcnta (isldr, isldr1, xsldr, ntdim, 2)
            go to 4
         endif
         nptr(m) = k
         read (nfile2) (isldr(j + k),isldr1(j + k),xsldr(j + k),j=1,ncsf)
!          write(*,*) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
         k = k + ncsf
!
      end do
!
!   Close (and hence release) the scratch file
!
!     IF (AVAIL) THEN
!
!        CLOSE (unit=NFILE,status="DELETE")
!      ENDIF
!
      nptr(m+1) = k
      ninteg = m
!
      return
!
  301 format(/,/,/,1x,i8,' MCT coefficients generated for rank ',i2,&
         ' and parity ',i2,/,/)
      return
!
      end subroutine mctin
