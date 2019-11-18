module setup_parameter
    use fileIO
    implicit none
    public

    ! ## declare user setting parameter ##
    integer, parameter :: nair = 5
    integer, parameter :: npml = 12
    integer, parameter :: thread_num = 4

    integer, parameter :: nx = 40 + (nair + npml) * 2
    integer, parameter :: ny = 66 + (nair + npml) * 2
    integer, parameter :: nz = 40 + (nair + npml) * 2
    integer, parameter :: nt = 30001
    integer, parameter :: check_time = 20
    integer, parameter :: check_interval = int(nt / check_time)
    
    double precision, parameter :: freq = 6.78d+6
    double precision, parameter :: dx = 0.002d0, dy = dx, dz = dx

    integer, parameter :: dimpml = 3
    double precision, parameter :: refpml = 1.0d-6


    ! ## declare general parameter ##
    double precision, parameter :: pi = 3.1415926535d0
    double precision, parameter :: hpi = pi / 2.0d0
    double precision, parameter :: cv = 2.99792458d+8
    double precision, parameter :: eps0 = 8.8541878128d-12
    double precision, parameter :: mu0 = pi * 4.0d-7
    double precision, parameter :: z0 = 120 * pi
    double precision, parameter :: omega = 2.0d0 * pi * freq
    double precision, parameter :: lambda = cv/freq
    double precision, parameter :: wavenum = 2.0d0 * pi / lambda
    complex(kind(0d0)), parameter :: im = (0.0d0, 1.0d0)

    integer, parameter :: cent_x = int(nx/2)
    integer, parameter :: cent_y = int(ny/2)
    integer, parameter :: cent_z = int(nz/2)

    double precision :: dt
    double precision :: hdt
    
    complex(kind(0d0)) :: cphase

    ! ## test model para ##
    integer, parameter :: mie_radius = int(0.04d0 / dx)
    integer, parameter :: mie_per = 2
    integer, parameter :: mie_distance = int(0.01d0 / dy)
    integer, parameter :: mie_dipole_len = int(16)

    integer, parameter :: feedx = cent_x
    integer, parameter :: feedy = cent_y + mie_radius + mie_distance
    integer, parameter :: feedz = cent_z

    integer, parameter :: ground_plane_height = npml + 1


    ! ## declare public variable ##
    integer :: step = 0
    integer :: ios = 1
    double precision :: t = 0
    double precision :: vfeed, ifeed
    double precision :: vfeed_sub


    ! ## declare EM-parameter array ##
    integer, parameter :: nmax_per = 5

    integer, dimension(nx, ny, nz) :: idper = 1

    integer, dimension(nx, ny, nz) & 
        :: idperx = 1, idpery = 1, idperz = 1 ! permitivity ID mapping on calc-space for each axis
    integer, dimension(nx, ny, nz) & 
        :: idpecx = 1, idpecy = 1, idpecz = 1 ! PEC ID mapping on calc-space for each axis

    double precision, dimension(nmax_PER) :: sigma = 0.0d0
    double precision, dimension(nmax_PER) :: eps = eps0
    double precision, dimension(nmax_PER) :: msigma = 0.0d0
    double precision, dimension(nmax_PER) :: mu = mu0
    double precision, dimension(nmax_PER) :: rho = 1.0d0


    !## declare EM-field array ##
    double precision, dimension(nx, ny, nz) :: ex, ey, ez
    double precision, dimension(nx, ny, nz) :: hx, hy, hz

    double precision, dimension(nx, ny, nz) :: ex_sub, ey_sub, ez_sub
    double precision, dimension(nx, ny, nz) :: hx_sub, hy_sub, hz_sub

    double precision, dimension(nx, ny, nz) :: examp, eyamp, ezamp
    double precision, dimension(nx, ny, nz) :: hxamp, hyamp, hzamp

    double precision, dimension(nx, ny, nz) :: sar
    double precision :: sar_ave_wb = 0.d0
    double precision :: mass_weight = 0.d0

    double precision, dimension(nx, ny, nz) :: exphase, eyphase, ezphase
    double precision, dimension(nx, ny, nz) :: hxphase, hyphase, hzphase

    complex(kind(0d0)), dimension(nx, ny, nz) :: einx, einy, einz
    complex(kind(0d0)), dimension(nx, ny, nz) :: einx_sub, einy_sub, einz_sub

    double precision, dimension(nx, ny, nz) :: etx, ety, etz
    double precision, dimension(nx, ny, nz) :: etx_sub, ety_sub, etz_sub

    !## declare maxwell's co-efficient array ##
    double precision, dimension(nx, ny, nz) :: cex, cey, cez
    double precision, dimension(nx, ny, nz) :: dex, dey, dez

    double precision, dimension(nx, ny, nz) :: chx, chy, chz
    double precision, dimension(nx, ny, nz) :: dhx, dhy, dhz

end module setup_parameter