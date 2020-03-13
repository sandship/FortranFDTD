program main
    implicit none
    call make_plane_wave

end program main


subroutine make_plane_wave
    use setup_parameter
    implicit none
    integer :: i, j, k

    open(101, file='incfield/planewave_xaxis_6_78MHz.dat', status='replace')
    do k = -nz/2, nz/2
        do j = -ny/2, ny/2
            do i = -nx/2, nx/2
                write(101, *) i * dx, j * dy, k * dz, & 
                              real(exp(im * (k * dz * wavenum))), &
                              imag(exp(im * (k * dz * wavenum))), &
                              0.d0, &
                              0.d0, &
                              0.d0, &
                              0.d0
            end do
        end do
    end do
    close(101)

end subroutine make_plane_wave