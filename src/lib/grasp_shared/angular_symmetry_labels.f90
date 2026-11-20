      module angular_symmetry_labels
      implicit none (type, external)
      character(len=4) :: jlbl(32), jlbr(32)
      character(len=4) :: jlbp(2)
      ! Labels indexed by 2J + 1, in left- and right-justified forms.
      data jlbl/"0   ", "1/2 ", "1   ", "3/2 ", "2   ", "5/2 ", "3   ", "7/2 ",&
                "4   ", "9/2 ", "5   ", "11/2", "6   ", "13/2", "7   ", "15/2",&
                "8   ", "17/2", "9   ", "19/2", "10  ", "21/2", "11  ", "23/2",&
                "12  ", "25/2", "13  ", "27/2", "14  ", "29/2", "15  ", "31/2"/
      data jlbr/"   0", " 1/2", "   1", " 3/2", "   2", " 5/2", "   3", " 7/2",&
                "   4", " 9/2", "   5", "11/2", "   6", "13/2", "   7", "15/2",&
                "   8", "17/2", "   9", "19/2", "  10", "21/2", "  11", "23/2",&
                "  12", "25/2", "  13", "27/2", "  14", "29/2", "  15", "31/2"/
      ! Negative and positive parity labels, in that order.
      data jlbp/ " -  ", " +  "/

      end module angular_symmetry_labels
