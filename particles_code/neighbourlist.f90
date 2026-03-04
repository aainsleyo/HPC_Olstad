subroutine ONETO2(id, x, y)
    implicit none
    integer, intent(in)  :: id
    integer, intent(out) :: x, y

    y = id / b
    x = mod(id, b)
end subroutine ONETO2

function TWOtoONE(x, y) result(id)
    implicit none
    integer, intent(in) :: x, y
    integer :: id

    id = y * b + x
end function TWOtoONE

function get_neighbour_ids(p) result (neighbours)
    integer, intent(in) :: p !> The sequential position
    integer, allocatable, dimension(:) :: neighbours !> The list of neighbours
    integer :: i,j,k
    integer :: x_cell, y_cell, x_test, y_test
    integer :: max_list(9)

    !> Don't save the values between calls
    k = 1
    max_list = -1

    call ONETO2(p,x_cell,y_cell)
    max_list(k) = p

    !print *, max_list

    !> Start by just getting all of the neighbours
    do i = -1,1
        do j = -1,1
            if(i == 0 .and. j == 0) cycle !> Skip the cell itself
            
            !> The "test" coordinate
            x_test = x_cell + i
            y_test = y_cell + j

            !> If the coordinates are real (within the expected range), add them to the list:
            if((x_test .ge. 0 .and. x_test .lt. b) .and. (y_test .ge. 0 .and. y_test .lt. b)) then
                !> Increment the number of correct coordinates
                k = k + 1
                !> Flip a zero to the new coordinate in the max list
                max_list(k) = TWOtoONE(x_test, y_test)
            end if
        end do
    end do

    !> Return a vector of all non-zero elements
!    neighbours = pack(max_list, max_list>= 0.0)
    neighbours = max_list ! entries -1 must be skipped as they correspond to non-existent neighbours at the edge or corner
end function get_neighbour_ids
    
subroutine buildNBL()
  integer :: sector
  do sector=0,b*b-1
     nbl(sector,:) = get_neighbour_ids(sector)
  end do
end subroutine buildNBL
