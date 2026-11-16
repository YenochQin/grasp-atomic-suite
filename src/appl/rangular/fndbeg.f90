!***********************************************************************
!                                                                      *
      subroutine fndbeg(jastrt, jbstrt, index, llistt, llistv)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use mcp_c
      use orb_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: jastrt
      integer, intent(inout) :: jbstrt
      integer, intent(inout) :: index
      integer, intent(out) :: llistt
      integer, dimension(:), pointer :: llistv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, ios, icread, irread, indx, nrec, i, label
      real(double) :: coeff
      character :: srtlab*8, mcplab*3
!-----------------------------------------------
!
!   Begin by examining file 30; this is the last to be updated by
!   SUBROUTINE MCP
!
!   Read and check the character part of the file header
!
      rewind (30)
      read (30) mcplab, srtlab
      if (srtlab == '  SORTED') then
         jastrt = ncf + 1
         jbstrt = ncf + 1
         llistt = 0
         llistv(:kmax) = 0
         go to 8
      endif
!
      read (30)
      read (30)
!
!   Read as many records as possible
!
    2 continue
      read (30, iostat=ios) icread, irread, indx
!
      if (ios == 0) then
!
!   No errors or end-of-file; keep reading
!
         jastrt = icread
         jbstrt = irread
         index = indx
         go to 2
!
      else
!
         if (jastrt==ncf .and. jbstrt==ncf) then
!
!   All coefficients have been generated; sorting may still
!   be necessary; force this option
!
            jastrt = ncf + 1
            jbstrt = ncf + 1
!
         else
!
!   Some coefficients remain to be generated; reposition all files
!   for augmentation of lists by SUBROUTINE MCP; update JBSTRT and,
!   if appropriate, JASTRT
!
            do k = 31, 32 + kmax
               rewind (k)
               nrec = 3
               do i = 1, nrec
                  read (k)
               end do
    4          continue
               read (k, iostat=ios) indx, label, coeff
               if (ios==0 .and. indx<=index) then
                  nrec = nrec + 1
                  go to 4
               else
                  rewind (k)
                  do i = 1, nrec
                     read (k)
                  end do
                  if (k > 31) then
                     llistv(k-32) = nrec - 3
                  else
                     llistt = nrec - 3
                  endif
               endif
            end do
!
!  Now, reposition the sms file. This file should contain the
!  same number of data records as file 33.

            rewind (20)
            do i = 1, llistv(1)
               read (20)
            end do

            jbstrt = jbstrt + 1
            if (jbstrt > ncf) then
               jastrt = jastrt + 1
               jbstrt = jastrt
            endif
!
         endif
!
      endif
!
    8 continue
      return
      end subroutine fndbeg
