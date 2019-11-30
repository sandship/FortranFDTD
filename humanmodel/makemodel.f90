program main
    implicit none
    integer, parameter :: radius = 20
    integer, parameter :: cent_x = int(radius / 2)
    integer, parameter :: cent_y = int(radius / 2)
    integer, parameter :: cent_z = int(radius / 2)

    integer :: i, j, k
    open(101, file='test_sphere.index', status='replace')

    do i = 0, radius
        do j = 0, radius
            do k = 0, radius
                if((k - cent_z) ** 2 + (j - cent_y) ** 2 + (i - cent_x) ** 2 <= (radius * 0.5d0) ** 2) then
                    write(101, *) i, j, k, 3
                end if
            end do
        end do
    end do

    close(101)
    
end program main