!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**        *****   ******  **   **  **    **   *****   ******        ***
!**       **   **  **      ***  **  ***  ***  **   **  **   **       ***
!**       **       **      ***  **  ** ** **  **       **   **       ***
!**       **  ***  ****    ** ****  ** ** **  **       ******        ***
!**       **   **  **      **  ***  **    **  **       **            ***
!**       **   **  **      **   **  **    **  **   **  **            ***
!**        *****   ******  **   **  **    **   *****   **            ***
!**                                                                  ***
!**   Program for generating the energy expression for H(DC). This   ***
!**   program is a derivative of GRASP2 (F. A. Parpia, I. P. Grant,  ***
!**   and C. F. Fischer, 1990).                                      ***
!**                                                                  ***
!**                         GRASP92 Version                          ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      PROGRAM GENMCPMPI
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
!   Updated by Charlotte F. Fischer
!                                                                      *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
      USE parameter_def,   ONLY:  NNNW
      USE iso_fortran_env, ONLY: int64
      USE rangular_mpi_input, ONLY: mpi_angular_input_options, parse_mpi_angular_options, &
           open_mpi_angular_profile, load_mpi_angular_header, load_mpi_angular_block, close_mpi_angular_input
      USE rangular_profile, ONLY: stage_clock, record_stage, close_stage_profile
      USE memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE mpi_C
      USE default_C
      Use hblock_c
      Use orb_C
      Use stat_c
      Use mcp_c
      Use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbg_I
      USE setmc_I
      USE setsum_I
      USE cslh_I
      USE strsum_I
      USE factt_I
      USE settmp_I
      USE mpiu_I, ONLY: startmpi2
      USE mpiu_I, ONLY: stopmpi2
      USE setdbgmpi_I, ONLY: setdbgmpi
      USE setmcpmpi_I, ONLY: setmcpmpi
      USE mcpmpi_I, ONLY: mcpmpi
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER NBLK0, NCOUNT1, NCORE, NB
      PARAMETER (NBLK0 = 50)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL ::  DEBUG, RESTRT, YES
      CHARACTER(LEN=128) :: STARTDIR, PERMDIR, TMPDIR
      CHARACTER(LEN=8), DIMENSION(NBLK0):: IDBLK
      CHARACTER(LEN=3) :: IDSTRING
      INTEGER :: LENPERM, LENTMP
      TYPE(mpi_angular_input_options) :: input_options
      INTEGER(int64) :: program_start_count, stage_start_count
      CHARACTER(len=16) :: input_format_name
!-----------------------------------------------
      CALL parse_mpi_angular_options(input_options)
      OPEN(UNIT=739,FILE='rangular.log',STATUS='UNKNOWN')
      startdir = ' '
      permdir = ' '
      tmpdir = ' '
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1, &
                           startdir, permdir, tmpdir, 'RANGULAR_MPI')
      WRITE (idstring, '(I3.3)') myid
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)
      CALL open_mpi_angular_profile(input_options)
      input_format_name = 'text'
      IF (input_options%use_descriptors) input_format_name = 'descriptors'
      program_start_count = stage_clock()
!=======================================================================
!  Get NDEF on node-0 and then send to all nodes
!=======================================================================
      IF (myid == 0) THEN
         WRITE (istde,*)
!         WRITE (istde,'(A)') ' Default settings?  (y/n) '
         WRITE (istde,'(A)') ' Full interaction?  (y/n) '
         YES = GETYN ()
         IF (YES) THEN
            NDEF = 0
            write(739,'(A)') 'y            ! Full interaction'
         ELSE
            NDEF = 1
            write(739,'(A)') 'n            ! Full interaction'
         ENDIF
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    setdbg - debug output control parameters
!    setmc - machine- and installation- dependent constants
!    setsum - open the summary file
!    cslh - load header of the csl file
!    setmcp - open and check the  .mcp  files
!    strsum - append a summary of the inputs to the  .sum  file
!    factt - table of logarithms of factorials setup
!=======================================================================

      CALL setdbgmpi (DEBUG, permdir(1:lenperm) // '/genmcp.dbg')
      CALL SETMC
      stage_start_count = stage_clock()
      CALL load_mpi_angular_header(input_options,permdir,ncore,nblk0,idblk)
      CALL record_stage('input_header_and_broadcast',0,stage_start_count,storage=input_format_name)
      stage_start_count = stage_clock()
      RESTRT = .FALSE.
!cjb  myid, nprocs = NOT args
!cjb  CALL setmcpmpi (myid, nprocs, ncore, idblk, 'mcp' // idstring)
      CALL setmcpmpi (ncore, idblk, 'mcp' // idstring)
      IF (NDEF/=0 .AND. MYID==0) CALL STRSUM
      CALL FACTT
      CALL record_stage('mcp_setup_including_interaction_input',0,stage_start_count)

!=======================================================================
!     For each block, generate and sort the data
!=======================================================================

      DO NB = 1, NBLOCK
         NCF = NCFBLK(NB)                        ! This ncf goes to common
         IF (MYID == 0) THEN
            WRITE (6, *)
            WRITE (6, *) 'Block ', NB, ',  ncf = ', NCF
         ENDIF
            !*** Load current CSL block. Memories de-allocated in mcp ***
         stage_start_count = stage_clock()
         CALL ALLOC (iqa, NNNW, NCF, 'IQA', 'GENMCP')
         CALL ALLOC (jqsa, NNNW,3,NCF, 'JQSA', 'GENMCP')
         CALL ALLOC (jcupa, NNNW, NCF, 'JCUPA', 'GENMCP')
         CALL record_stage('csf_allocation',nb,stage_start_count,int(ncf,int64))
!
         stage_start_count = stage_clock()
         CALL load_mpi_angular_block(input_options,ncore,nb)
         CALL record_stage('csf_load_and_broadcast',nb,stage_start_count,int(ncf,int64),storage=input_format_name)
         !*** Open tmp.xx files for block nb ***
         CALL SETTMP (NB, KMAX, 'tmp' // idstring)
         !*** Generation of MCP coefficients ***
!cjb  myid, nprocs = NOT args
!cjb     CALL MCPmpi (NB, RESTRT, MYID, NPROCS, 'mcp' // idstring)
         CALL MCPmpi (NB, RESTRT, 'mcp' // idstring)
      END DO
      CALL close_mpi_angular_input(input_options)
      CLOSE(24)                                  ! Summary file
      CLOSE(739)                                 ! rangular.log
      IF (DEBUG) CLOSE(99)                       ! Debug file
      CALL record_stage('program_before_finalize_including_interaction_input',0,program_start_count)
      CALL close_stage_profile()
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      CALL stopmpi2 (myid,nprocs,host,lenhost,ncount1,'RANGULAR_MPI')
      STOP
      END PROGRAM GENMCPMPI
