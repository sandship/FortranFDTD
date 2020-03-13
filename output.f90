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

    ! at first
    subroutine output_at_first
        implicit none
        integer :: i, j, k
        character(len=9) :: stepchar
        write(stepchar, '(i9.9)') step

        open(101, file='./output/load_planewave.dat', status='replace')
        i = cent_x
        do j = 1, ny
            do k = 1, nz
                write(101, *) j, k, real(einx(i, j, k)), imag(einx(i, j, k))
            end do
        end do
        close(101)

        open(101, file='./output/CE_XYplane.dat', status='replace')
            k = cent_z
            do i = 1, nx
                do j = 1, ny
                    write(101, *) i, j, cex(i, j, k), cey(i, j, k), cez(i, j, k) &
                                      , dex(i, j, k), dey(i, j, k), dez(i, j, k)
                end do
            end do
        close(101)
    
        open(101, file='./output/CE_YZplane.dat', status='replace')
            i = cent_x
            do k = 1, nz
                do j = 1, ny
                    write(101, *) j, k, cex(i, j, k), cey(i, j, k), cez(i, j, k) &
                                      , dex(i, j, k), dey(i, j, k), dez(i, j, k)
                end do
            end do
        close(101)
    
        open(101, file='./output/CE_XZplane.dat', status='replace')
            j = cent_y
            do k = 1, nz
                do i = 1, nx
                    write(101, *) i, k, cex(i, j, k), cey(i, j, k), cez(i, j, k) &
                                      , dex(i, j, k), dey(i, j, k), dez(i, j, k)
                end do
            end do
        close(101)
    end subroutine output_at_first


    ! at the end
    subroutine output_at_end
        implicit none
        integer :: n
        character(len=9) :: stepchar
        write(stepchar, '(i9.9)') step

    end subroutine output_at_end


    ! at middle
    subroutine output_at_checkpoint
        implicit none
        integer :: i, j, k
        character(len=9) :: stepchar
        write(stepchar, '(i9.9)') step

        open(101, file='./output/Eamp_center_timeline.dat', position='append')
            i = cent_x
            j = cent_y
            k = cent_z
            write(101, *) step, sqrt(examp(i, j, k)**2 + eyamp(i, j, k)**2 + ezamp(i, j, k)**2) &
                              , sqrt(etx_sub(i, j, k)**2 + ety_sub(i, j, k)**2 + etz_sub(i, j, k)**2) &
                              , sqrt(etx(i, j, k)**2 + ety(i, j, k)**2 + etz(i, j, k)**2)
        close(101)

        open(101, file='./output/SARwb_timeline.dat', position='append')
            i = cent_x
            j = cent_y
            k = cent_z
            write(101, *) step, sar_ave_wb
        close(101)


        open(101, file='./output/SAR_YZplane'// stepchar //'.dat', status='replace')
            i = cent_x
            do k = npml, nz - npml
                do j = npml, ny - npml
                    write(101, *) j, k, sar(i, j, k)
                end do
            end do
        close(101)

        open(101, file='./output/SAR_XYplane'// stepchar //'.dat', status='replace')
            k = cent_z
            do j = npml, ny - npml
                do i = npml, nx - npml
                    write(101, *) i, j, sar(i, j, k)
                end do
            end do
        close(101)

        open(101, file='./output/SAR_XZplane'// stepchar //'.dat', status='replace')
            j = cent_y
            do k = npml, nz - npml
                do i = npml, nx - npml
                    write(101, *) i, k, sar(i, j, k)
                end do
            end do
        close(101)

        open(101, file='./output/Eamp_YZplane'// stepchar //'.dat', status='replace')
            i = cent_x
            do k = npml, nz - npml
                do j = npml, ny - npml
                    write(101, *) j, k, examp(i, j, k), eyamp(i, j, k), ezamp(i, j, k) &
                                      , sqrt(examp(i, j, k)**2 + eyamp(i, j, k)**2 + ezamp(i, j, k)**2)
                end do
            end do
        close(101)
    
        open(101, file='./output/Eamp_XYplane'// stepchar //'.dat', status='replace')
            k = cent_z
            do j = npml, ny - npml
                do i = npml, nx - npml
                    write(101, *) i, j, examp(i, j, k), eyamp(i, j, k), ezamp(i, j, k) &
                                      , sqrt(examp(i, j, k)**2 + eyamp(i, j, k)**2 + ezamp(i, j, k)**2)
                end do
            end do
        close(101)
    
        open(101, file='./output/Eamp_XZplane'// stepchar //'.dat', status='replace')
            j = cent_y
            do k = npml, nz - npml
                do i = npml, nx - npml
                    write(101, *) i, k, examp(i, j, k), eyamp(i, j, k), ezamp(i, j, k) &
                                      , sqrt(examp(i, j, k)**2 + eyamp(i, j, k)**2 + ezamp(i, j, k)**2)
                end do
            end do
        close(101)
    
    end subroutine output_at_checkpoint

end module output