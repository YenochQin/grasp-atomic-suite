!***********************************************************************
      subroutine lodcslmpigg (nfile, ncore, jblock)

! An MPI container of lodcsh2 which loads CSL list of the current block
! into memory. It forwards the call together with the same set of
! parameters to lodcsh2 and then broadcasts the results to all nodes.
!
! Note: Memories have been allocated/deallocated each block outside.
! This subroutine calls lodcsh2 on node-0 to generate the data for the
! block; and then broadcasts to all other nodes. A new MPI data type
! of 4 byte-long is created to handle 64-bit machines whose MPI
! implementation does not support 4-byte integers. If jblock=-119,
! then ALL blocks will be loaded instead of just one. This is
! implemented in lodcsh2.
!
! Currently used by rcimpivu, mcpmpi, rscfmpivu
!
! Xinghong He 98-08-06
!
!***********************************************************************
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,   only: max_orbitals
      use orb_c,           only: csf_count, csf_subshell_occupations
      use syma_c,          only: jpgg, nblk0
      use mpi_c
      use lodcsh2gg_i, only: lodcsh2gg
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer       :: nfile, ncore, jblock
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mpix_int1
!-----------------------------------------------------------------------
!
      if (mpi_rank == 0) then
!GG         CALL lodcsh2 ((nfile), (ncore), (jblock))
         call lodcsh2gg ((nfile), (ncore), (jblock))
      end if

! Construct mpi data type for Integer*1 and then broadcast.
!cjb mpix_bytes
!cjb if your compiler or your MPI does not accept type MPIX_INT1
!cjb try standard MPI type MPI_INTEGER1
      call mpix_bytes (1, mpix_int1, mpi_error_code)
      call mpi_bcast (csf_subshell_occupations(:,:),max_orbitals*csf_count,mpix_int1,0,mpi_comm_world,mpi_error_code)
!cjb  CALL MPI_Bcast (IQA(:,:),NNNW*NCF,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (jpgg(:),nblk0,mpi_integer,0,mpi_comm_world,mpi_error_code)
!GG      CALL MPI_Bcast                                                   &
!GG               (JQSA(:,:,:),3*NNNW*NCF,MPIX_INT1,0,MPI_COMM_WORLD,ierr)
!GG      CALL MPI_Bcast                                                   &
!GG               (JCUPA(:,:),   NNNW*NCF,MPIX_INT1,0,MPI_COMM_WORLD,ierr)
      return
      end subroutine lodcslmpigg
