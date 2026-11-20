!***********************************************************************
!                                                                      *
      subroutine spodmv(n, m, b, c)
!                                                                      *
!   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
!   lower triangle of the  (NxN)  matrix  A  is read from the disk     *
!                                                                      *
!   This is an adaptation of  Andreas Stathopulos'  routine  SPSBMV,   *
!   and is specific to GRASP2 derivatives.                             *
!                                                                      *
!   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
!               [SPBLAS]: DAXPYI/SAXPYI, DDOTI/SDOTI.                  *
!                                                                      *
!   F A Parpia and A Stathopoulos         Last revision: 13 Oct 1992   *
!   Block Version by Xinghong He          Last revision: 18 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use where_c, only: imcdf
      use rci_matrix_file_position, only: nposition
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use real_vector_initialization, only: dinit
      use file_positioning, only: posfile
      use sparse_dot_and_scatter, only: dmerge
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(real64), intent(in) :: b(n, m)
      real(real64), intent(inout) :: c(n, m)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: irow(n)
      integer :: myid, nprocs, ncf, ncfdum, iccutdum, myiddum, nprocsdum, icol&
         , nelc, ir, iv
      real(real64) :: emt(n)
      real(real64) :: elsto, diag, dl
!-----------------------------------------------
!
!     !...nposition+1 is the current position of the .res file
!     !...It is set in matrix and used in maneig, spodmv
!
!-----------------------------------------------------------------------
      write (6, *) "Calling spodmv..."
      myid = 0
      nprocs = 1
      ncf = n

!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON

      call dinit (n*m, 0.0_real64, c, 1)

      !...moved from maneig before "CALL GDVD (SPODMV..."
      call posfile (0, imcdf, nposition)

      read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
      if (ncf/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) stop &
         "spodmv: ncf read wrong"


      do icol = myid + 1, n, nprocs
         read (imcdf) nelc, elsto, (emt(ir),ir=1,nelc), (irow(ir),ir=1,nelc)
         do iv = 1, m
            diag = c(icol,iv) + (emt(nelc)-hamiltonian_energy_shift)*b(icol,iv)
            call dmerge (nelc - 1, b(1,iv), c(1,iv), irow(1), emt(1), b(icol,iv&
               ), dl)
            c(icol,iv) = diag + dl
         end do
      end do


      return
      end subroutine spodmv
