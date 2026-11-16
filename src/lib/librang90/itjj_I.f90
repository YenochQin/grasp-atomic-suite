      module itjj_i
      interface
!
      integer function itjj(ik,id,kg,bk,bd,ibt,bt,kg1,itp,itg,iq)
      use vast_kind_param, only:  double
      integer,      intent(in)                :: kg, kg1, iq
      integer,      intent(out)               :: itp, itg
      integer,      intent(in),  dimension(7) :: ik, id
      integer,      intent(out), dimension(7) :: ibt
      real(double), intent(in),  dimension(3) :: bk, bd
      real(double), intent(out), dimension(3) :: bt
      end function
      end interface
      end module
