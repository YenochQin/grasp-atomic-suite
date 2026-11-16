module itjpo_i
    implicit none (type, external)
    private
    public :: itjpo
contains
integer function itjpo(index)
    implicit none (type, external)
    integer,intent(in) :: index
    if(index<1.or.index>3) error stop 'Unexpected level index'
    itjpo=2*index-1
end function
end module itjpo_i
