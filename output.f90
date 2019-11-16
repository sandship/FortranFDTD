module output
    use fileIO
    use setup_parameter
    use load_model
    use load_field
    use field_update
    use calc_amp
    use calc_scatter
    implicit none


    contains

    subroutine output_at_first
        implicit none
        integer :: i, j, k

        open(101, file='./output/CE_YZplane.dat', status='replace')
        i = cent_x
        do k = npml, nz - npml
            do j = npml, ny - npml
                write(101, *) j, k, cex(i, j, k), cey(i, j, k), cez(i, j, k) &
                                  , dex(i, j, k), dey(i, j, k), dez(i, j, k)
            end do
        end do
        close(101)
    
        open(101, file='./output/CE_XZplane.dat', status='replace')
        j = cent_y
        do k = npml, nz - npml
            do i = npml, nx - npml
                write(101, *) i, k, cex(i, j, k), cey(i, j, k), cez(i, j, k) &
                                  , dex(i, j, k), dey(i, j, k), dez(i, j, k)
            end do
        end do
        close(101)
    end subroutine output_at_first


    subroutine output_at_end
        implicit none
        integer :: i, j, k

        open(101, file='./output/Eamp_Yline.dat', status='replace')
        i = cent_x
        k = cent_z
        do j = npml, ny - npml
            write(101, *) j, examp(i, j, k), eyamp(i, j, k), ezamp(i, j, k) &
                           , sqrt(examp(i, j, k)**2 + eyamp(i, j, k)**2 + ezamp(i, j, k)**2)
        end do
        close(101)
    
    end subroutine output_at_end


    subroutine output_at_checkpoint
        implicit none
        integer :: i, j, k
        character(len=8) :: stepchar
        write(stepchar, '(I8)') step

        open(101, file='./output/'// stepchar //'_Eamp_YZplane.dat', status='replace')
        i = cent_x
        do k = npml, nz - npml
            do j = npml, ny - npml
                write(101, *) j, k, examp(i, j, k), eyamp(i, j, k), ezamp(i, j, k) &
                , sqrt(examp(i, j, k)**2 + eyamp(i, j, k)**2 + ezamp(i, j, k)**2)
            end do
        end do
        close(101)

    end subroutine output_at_checkpoint

end module output