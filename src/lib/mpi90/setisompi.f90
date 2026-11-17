module setisompi_i
    implicit none (type, external)
    private
    public :: setisompi
contains
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mpi_c, only: double, mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_integer, mpi_rank
      use def_c, only: double, emn, nuclear_charge
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use nsmdat_c, only: sqn, dmomnm, qmomb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setiso_i, only: setiso
      use orthsc_i, only: orthsc
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
      if (mpi_rank == 0) call setiso (isofile)

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
      end subroutine setisompi
end module setisompi_i
