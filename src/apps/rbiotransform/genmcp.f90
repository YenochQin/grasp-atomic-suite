!***********************************************************************
!                                                                      *
      subroutine genmcp(name, ic, ntestg, inpci)
!                                                                      *
!   Entry routine for GENMCP. Controls the computation of the          *
!   one-particle coupling coefficients.                                *
!                                                                      *
!   Written by Per Jonsson                                 June 1996   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcslb_i, only: setcslb
      use log_factorial_table, only: factt
      use mcpout_i, only: mcpout
      use mcpin_i, only: mcpin
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ic
      integer  :: ntestg
      integer  :: inpci
      character  :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: unused_core_subshell_count
!-----------------------------------------------
!
      write (6, *)
      write (6, *) "GENMCP: Execution begins for ", name
!
!   Open, check, load data from, and close the  csl  file
!
      call setcslb (name, unused_core_subshell_count,3)
!
!   Set up the table of logarithms of factorials for use by
!   angular modules
!
      call factt
!
!   Proceed with the generation of MCP coefficients
!
      call mcpout (name, ic, ntestg, inpci)
      call mcpin (name, ic, ntestg, inpci)
!
!   Print completion message
!
      write (6, *)
!
      return
      end subroutine genmcp
