!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**           **   **  *******   *****    *****    *****             ***
!**           **   **  **       **   **  **   **  **   **            ***
!**           **   **  **       **       **   **      **             ***
!**           *******  ****      *****    *****      **              ***
!**           **   **  **            **      **     **               ***
!**           **   **  **       **   **     **     **                ***
!**           **   **  **        *****    **      *******            ***
!**                                                                  ***
!**            Relativistic Hyperfine Structure Program              ***
!**                         GRASP92 Version                          ***
!**                         Dynamic Storage                          ***
!**                                                                  ***
!**   ============================================================   ***
!**   Copyright (c) 1995 by P Jonsson, F A Parpia, and C F Fischer   ***
!**   ============================================================   ***
!**   All rights reserved.  No part of this software or its accom-   ***
!**   panying documentation may be reproduced, stored in a retrie-   ***
!**   val system,  or transmitted,  in any form or  by any  means,   ***
!**   electronic, mechanical,  photocopying,  recording, or other-   ***
!**   wise, without the prior written permission of the authors.     ***
!**                                                                  ***
!**                           Disclaimer                             ***
!**                           ==========                             ***
!**   The  authors make  no warranties,  express or implied,  that   ***
!**   this software or its  accompanying documentation are free of   ***
!**   error or that  they will meet your requirements for any par-   ***
!**   ticular application.  Neither this software nor its accompa-   ***
!**   nying documentation should be relied  upon for solving prob-   ***
!**   lems if an incorrect solution could result in injury or loss   ***
!**   of, or damage to, property. If you use this software or  its   ***
!**   accompanying documentation,  you do so entirely  at your own   ***
!**   risk;  the authors disclaim all liability for direct or con-   ***
!**   sequential damage.                                             ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
      program hfs92
!                                                                      *
!   Entry routine for HFS92. Controls the entire computation.          *
!                                                                      *
!   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
!               [HFS92]: CHKPLT, GETHFD, HFS, SETDBG, SETSUM,          *
!                        STRSUM.                                       *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                         Last revision:    Nov 2017   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use default_c, only: double, inpci, ndef, ndump
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use setdbg_i, only: setdbg
      use setmc_i, only: setmc
      use setcon_i, only: setcon
      use setsum_i, only: setsum
      use setcsla_i, only: setcsla
      use gethfd_i, only: gethfd
      use getmixblock_i, only: getmixblock
      use strsum_i, only: strsum
      use factt_i, only: factt
      use hfsgg_i, only: hfsgg
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: k, nci, ncore_not_used
      logical   :: yes
      character :: name*24
!-----------------------------------------------
!

      write (stderr_unit, *)
      write (stderr_unit, *) "RHFS"
      write (stderr_unit, *) "This is the hyperfine structure program"
      write (stderr_unit, *) "Input files:  isodata, name.c, name.(c)m, name.w"
      write (stderr_unit, *) "Output files: name.(c)h, name.(c)hoffd"

      write (stderr_unit, *)
      write (stderr_unit, *) "Default settings?"
      yes = getyn()
      write (stderr_unit, *)
      if (yes) then
         ndef = 0
      else
         ndef = 1
      end if

   10 continue
      write (stderr_unit, *) "Name of state"
      read (*, "(A)") name
      k = index(name," ")
      if (k == 1) then
         write (stderr_unit, *) "Names may not start with a blank"
         goto 10
      end if
      write (stderr_unit, *)
      write (stderr_unit, *) "Mixing coefficients from a CI calc.?"
      yes = getyn()
      if (yes) then
         nci = 0
      else
         nci = 1
      end if
!
!   Check compatibility of plant substitutions
!
!GG      CALL CHKPLT
!
!   Determine if there is to be any debug printout; this will be
!   made on the  .dbg  file
!
      call setdbg
!
!   Perform machine- and installation-dependent setup
!
      call setmc
!
!   Set up the physical constants
!
      call setcon
!
!   Open the  .sum  file
!
      call setsum (name, nci)
!
!   Open, check, load data from, and close, the  .csl  file
!
      call setcsla (name, ncore_not_used)
!
!   Get the remaining information
!
      call gethfd (name)
!
!   Get the eigenvectors
!
!GG      WRITE (ISTDE, *) 'Block format?'
!GG      YES = GETYN()
!GG      WRITE (ISTDE, *)
!GG      IF (YES) THEN
         call getmixblock (name, nci)
!GG      ELSE
!GG         IF (NCI == 0) THEN
!GG            CALL GETMIXC (NAME)
!GG         ELSE
!GG            CALL GETMIXA (NAME)
!GG         ENDIF
!GG      ENDIF
!
!   Append a summary of the inputs to the  .sum  file
!
      call strsum
!
!   Set up the table of logarithms of factorials
!
      call factt
!
!   Proceed with the HFS calculation
!
      call hfsgg
!
!   Print completion message
!
      write (stderr_unit, *)
      write (stderr_unit, *) "RHFS: Execution complete."
!
      stop
      end program hfs92
