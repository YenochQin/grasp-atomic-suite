!***********************************************************************
      subroutine setisompi (isofile)

! An MPI container for setiso which opens, checks, and loads isofile
! to get isotope data. Data loaded are:
!     EMN,Z,/NPAR/,/NSMDAT/
! where /.../ means whole common block.
!
! Xinghong He 98-08-06
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mpi_c
      use def_c
      use npar_c
      use nsmdat_c, only: sqn, dmomnm, qmomb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setiso_i
      use orthsc_i
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = *) :: isofile
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, k, nwin, ios, npy, naky, my
!     integer :: ierr
      real(real64) :: con, fkk, ey, pzy, dnorm
      real(real64), pointer :: pa(:), qa(:), ra(:)
!-----------------------------------------------
      if (mpi_rank .eq. 0) call setiso (isofile)

      call mpi_bcast (nuclear_charge,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (emn,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (nuclear_model_parameters,2,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (nuclear_model_parameter_count,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (sqn,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (dmomnm,1,mpi_double_precision,0,mpi_comm_world, &
        mpi_error_code)
      call mpi_bcast (qmomb,1,mpi_double_precision,0,mpi_comm_world, &
        mpi_error_code)

      return
      end
