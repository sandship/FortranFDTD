module setup_parameter
    use fileIO
    implicit none
    public

    ! ## declare user setting parameter ##
    character(len=*), parameter :: human_model = "human_model/test_sphere.index"
    character(len=*), parameter :: inc_efield = "inc_field/planewave_xaxis_6_78MHz.dat"
    character(len=*), parameter :: tissue_param = "tissue_param/6_78MHZ_Tissue_param.csv"
    character(len=*), parameter :: antenna_name = "antenna_model/dipole"

    integer, parameter :: nair = 15
    integer, parameter :: npml = 12

    integer, parameter :: nx = 20 + (nair + npml) * 2
    integer, parameter :: ny = 20 + (nair + npml) * 2
    integer, parameter :: nz = 20 + (nair + npml) * 2
    integer, parameter :: nt = 2000
    integer, parameter :: check_time = 10
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
    double precision, parameter :: z0 = 120.0d0 * pi
    double precision, parameter :: omega = 2.0d0 * pi * freq
    double precision, parameter :: lambda = cv / freq
    double precision, parameter :: wn = 2.0d0 * pi / lambda !wavenumber
    complex(kind(0d0)), parameter :: im = (0.0d0, 1.0d0)

    integer, parameter :: cent_x = int(nx/2)
    integer, parameter :: cent_y = int(ny/2)
    integer, parameter :: cent_z = int(nz/2)
    integer, parameter :: margin = (nair + npml)

    double precision :: dt
    double precision :: hdt
    
    complex(kind(0d0)) :: cphase

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
    double precision :: sar_ave_wb = 0.d0


    !## declare maxwell's co-efficient array ##
    double precision, dimension(nx, ny, nz) :: cex, cey, cez
    double precision, dimension(nx, ny, nz) :: dex, dey, dez

    double precision, dimension(nx, ny, nz) :: chx, chy, chz
    double precision, dimension(nx, ny, nz) :: dhx, dhy, dhz


    !## declare EM-field array ##
    double precision, dimension(nx, ny, nz) :: ex, ey, ez
    double precision, dimension(nx, ny, nz) :: hx, hy, hz

    double precision, dimension(nx, ny, nz) :: ex_sub, ey_sub, ez_sub
    double precision, dimension(nx, ny, nz) :: hx_sub, hy_sub, hz_sub

    double precision, dimension(nx, ny, nz) :: examp, eyamp, ezamp
    double precision, dimension(nx, ny, nz) :: hxamp, hyamp, hzamp

    double precision, dimension(nx, ny, nz) :: eamp
    double precision, dimension(nx, ny, nz) :: sar

    double precision, dimension(nx, ny, nz) :: exphase, eyphase, ezphase
    double precision, dimension(nx, ny, nz) :: hxphase, hyphase, hzphase

    complex(kind(0d0)), dimension(nx, ny, nz) :: einx, einy, einz
    complex(kind(0d0)), dimension(nx, ny, nz) :: einx_sub, einy_sub, einz_sub

    double precision, dimension(nx, ny, nz) :: etx, ety, etz
    double precision, dimension(nx, ny, nz) :: etx_sub, ety_sub, etz_sub

    !## declare scatter E-field param

    ! pp(np0)			[r]   構成点データ
    ! tl(2,nla)			[i]   三角を構成する点の番号データ
    ! erad(nla)			[r]   線の半径
    ! lrp(nla)			[i]   ＋の線リスト
    ! lrm(nla)			[i]   －の線リスト
    ! edp(nedp)			[i]   エッジを構成する点の番号
    ! freen(2,nedp)		[i]   エッジから浮遊する2点の番号
    ! inner(6,nta)		[i]   寄与する線の番号
    ! length(nla)		[r]   線の長さ
    ! center(3,nla)		[r]   線の中心点
    ! nlp(nla)			[i]	  ひとつの線が寄与するエッジの数
    ! inds(nant)		[i]   給電点のインデックス
    
    complex(kind(0d0)), dimension(nx, ny, nz) :: cur_x, cur_y, cur_z
    double precision, dimension(nx, ny, nz) :: xi, yi, zi

    integer :: vertex_num, seg_num, edge_num

    integer, allocatable :: tl(:, :), lrp(:), lrm(:), edp(:)
    integer, allocatable :: freen(:, :), seg_inner(:, :), nlp(:)
    integer, allocatable :: inds(:), ipiv(:), iw1(:)

    double precision, allocatable :: pp(:, :), erad(:), seg_length(:)
    double precision, allocatable :: seg_center(:, :)
    
    complex(kind(0d0)), allocatable :: escatter_x(:)
    complex(kind(0d0)), allocatable :: escatter_y(:)
    complex(kind(0d0)), allocatable :: escatter_z(:)

    complex(kind(0d0)), allocatable :: return_v(:)

end module setup_parameter