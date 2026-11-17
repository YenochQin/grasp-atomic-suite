!***********************************************************************
      subroutine outsdampi(lprint,nnonz, ncf)
!***********************************************************************
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use mpi_c, only: mpi_comm_world, mpi_error_code, mpi_integer, mpi_rank, mpi_reduce, mpi_sum
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nnonz
      integer, intent(in) :: ncf
      logical, intent(in) :: lprint
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=9), parameter :: myname ="outsdampi"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nnonz_a
!-----------------------------------------------
!
      call mpi_reduce (nnonz, nnonz_a, 1, mpi_integer, mpi_sum, 0,  &
               mpi_comm_world, mpi_error_code)
      if (mpi_rank == 0) &
      write(6, *)" ... complete; density of non-zero elements of H(DC):", &
         nnonz_a, "/", (ncf*(ncf + 1))/2
!
!   Debug printout
!
      if (lprint) then
        if (mpi_rank == 0) then
          write (99, *)
          write (99, *) "From ", myname, " :"
          write (99, 301) nnonz_a
          write (6, *) "This part not finished. See ", myname
          write (99, *) "This part not finished. See ", myname
        end if
      end if
  301 format(" Number of nonzero elements in H(DC): ",1i4)
  302 format(" Column ",1i2,", row ",1i2,", sparse matrix index ",1i4)
      return
      end subroutine outsdampi
