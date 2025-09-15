module m_pdsim
  use iso_fortran_env, only: error_unit, int64
  use m_interp_unstructured
  use m_particle_core
  use m_units_constants
  use m_config
  use m_cross_sec
  use m_gas
  use m_table_data
  use m_lookup_table

  implicit none
  private

  integer, parameter, public :: dp = kind(0.0d0)
  real(dp), parameter, public :: undefined_real = -1e100_dp
  character(len=*), parameter, public :: undefined_str = "NOT_SPECIFIED"

  ! Convert V/m to Townsend
  real(dp), parameter, public :: SI_to_Townsend = 1e21_dp

  ! Convert Townsend to V/m
  real(dp), parameter, public :: Townsend_to_SI = 1e-21_dp

  ! Ideal gas number density at 300 K and 1 bar
  real(dp), parameter :: N0_300K_1bar = 1e5_dp / (UC_boltzmann_const * 300)

  integer, parameter :: pdsim_coord_2d = 1
  integer, parameter :: pdsim_coord_3d = 2
  integer, parameter :: pdsim_coord_axi = 3

  !> Number of columns in transport data lookup table
  integer, parameter, public :: pdsim_ncols = 3

  !> Column in lookup table with ionization coefficient
  integer, parameter, public :: pdsim_col_alpha = 1

  !> Column in lookup table with attachment coefficient
  integer, parameter, public :: pdsim_col_eta = 2

  !> Column in lookup table with electron mobility coefficient
  integer, parameter, public :: pdsim_col_mu = 3

  real(dp), parameter :: material_threshold = 1e-8_dp

  !> Point data: K integral
  integer, public :: i_k_integral

  !> Point data: Kendall w-function
  integer, public :: i_w

  !> Point data: Kendall P(n_electron=0)
  integer, public :: i_p0

  !> Point data: Ionization coefficient
  integer, public :: i_alpha

  !> Point data: Attachment coefficient
  integer, public :: i_eta

  !> Point data: Effective ionization coefficient
  integer, public :: i_alpha_eff

  !> Point data: K^*, which is like K but for number of ionizations
  integer, public :: i_kstar

  !> Point data: Probability of no further ionization by first electron
  integer, public :: i_p_m1

  !> Point data: Travel time of avalanche
  integer, public :: i_avalanche_time

  !> Point data: Coordinates of avalanche arrival
  integer, public :: i_x1, i_x2, i_x3

  !> Point data: Travel time of ions
  integer, public :: i_ion_time

  !> Point data: Coordinates of positive ion arrival
  integer, public :: i_ion_x1, i_ion_x2, i_ion_x3

  !> Point data: Secondary emission coefficient for positive ions
  integer, public :: i_ion_gamma

  !> Point data: Estimated inception time
  integer, public :: i_inception_time

  !> Point data: Estimated inception probability
  integer, public :: i_inception_prob

  !> Point data: Estimated variance of inception probability
  integer, public :: i_inception_pvar

  !> Base name for output
  character(len=200), public, protected :: pdsim_output_name

  !> How verbose the simulations are
  integer, public, protected :: pdsim_verbosity

  !> How much output is written (0 means none)
  integer, public, protected :: pdsim_output_level

  !> Size for lookup tables used for interpolation
  integer, public, protected :: pdsim_table_size

  !> Positive ion mobility (m^2/Vs)
  real(dp), public, protected :: pdsim_ion_mobility

  !> Secondary emission coefficient for ions reaching a boundary. The first
  !> element is for reaching a domain boundary, the other elements are for
  !> the non-gas materials in increasing order.
  real(dp), public, protected, allocatable :: pdsim_ion_gamma(:)

  !> Whether there is a positive ion SEE coefficient
  logical, public, protected :: pdsim_ion_see_enabled

  !> Secondary emission coefficient for photons reaching a boundary. The first
  !> element is for reaching a domain boundary, the other elements are for
  !> the non-gas materials in increasing order.
  real(dp), public, protected, allocatable :: pdsim_photoemission_gamma(:)

  !> Whether there is a positive photon SEE coefficient on a boundary
  logical, public, protected :: pdsim_photoemission_enabled

  type(iu_grid_t), public    :: pdsim_ug
  integer, public, protected :: pdsim_pdata_field(3)
  integer, public, protected :: pdsim_icdata_material
  integer, public, parameter :: pdsim_gas_material_value = 0
  integer, public, protected :: pdsim_num_materials
  integer, public, protected :: pdsim_coord_system
  integer, public, protected :: pdsim_ndim
  logical, public, protected :: pdsim_axisymmetric

  type(LT_t), public, protected :: pdsim_tdtbl

  ! Public methods
  public :: pdsim_create_config
  public :: pdsim_initialize
  public :: pdsim_convert_r
  public :: pdsim_pointdata_average

contains

  !> Create configuration for pdsim module
  subroutine pdsim_create_config(cfg)
    type(CFG_t), intent(inout) :: cfg
    real(dp)                   :: dummy_real(0)

    call CFG_add(cfg, "input%mesh", undefined_str, &
         "Input mesh file in (.binda format)", required=.true.)
    call CFG_add(cfg, "input%axisymmetric", .false., &
         "Whether the mesh is axisymmetric (only for 2D)")
    call CFG_add(cfg, "input%coordinate_scale_factor", 1.0_dp, &
         "Scale factor for coordinates")
    call CFG_add(cfg, "input%field_scale_factor", 1.0_dp, &
         "Scale the electric field with this factor")
    call CFG_add(cfg, "input%field_component_names", [undefined_str], &
         "Names of the electric field components stored as point data", &
         dynamic_size=.true., required=.true.)
    call CFG_add(cfg, "input%material_name", undefined_str, &
         "Variable describing the material (point or cell data), or NONE", &
         required=.true.)
    call CFG_add(cfg, "input%lookup_table_size", 1000, &
         "Size to use for cross section lookup table")
    call CFG_add(cfg, "input%transport_data_file", undefined_str, &
         "File with electron transport and reaction data", required=.true.)
    call CFG_add(cfg, "input%input_interpolation", "linear", &
         "Input interpolation method (linear, cubic_spline)")
    call CFG_add(cfg, "input%interpolation_xspacing", "linear", &
         "Spacing used in lookup table (linear, quadratic, cubic)")
    call CFG_add(cfg, "input%three_body_first_species", 'O2', &
         "Primary species for three-body attachment")
    ! TODO: include more species here, like N2
    call CFG_add(cfg, "input%three_body_second_species", ['O2 ', 'H2O'], &
         "Third bodies for three-body attachment", &
         dynamic_size=.true.)
    call CFG_add(cfg, "input%three_body_efficiencies", [1.0_dp, 6.0_dp], &
         "Efficiencies of third bodies for three-body attachment", &
         dynamic_size=.true.)
    call CFG_add(cfg, "input%ion_mobility", 2e-4_dp * N0_300K_1bar, &
         "Reduced positive ion mobility [1/(m V s)]")
    call CFG_add(cfg, "input%ion_gamma", dummy_real, &
         "Secondary emission coefficient for ions reaching a boundary", &
         dynamic_size=.true.)
    call CFG_add(cfg, "input%photoemission_gamma", dummy_real, &
         "Secondary emission coefficient for photons reaching a boundary", &
         dynamic_size=.true.)

    call CFG_add(cfg, "gas%temperature", 300.0_dp, "Gas temperature (K)")
    call CFG_add(cfg, "gas%pressure", 1.0_dp, "Gas pressure (bar)")

    call CFG_add(cfg, "output%name", "output/pdsim", &
         "The base name for output files")
    call CFG_add(cfg, "output%verbosity", 1, &
         "How verbose the simulations are (higher means more output)")
    call CFG_add(cfg, "output%level", 2, &
         "How much output is written (0: none, 1: info, 2: all)")
    call CFG_add(cfg, "output%particles", .true., &
         "Output the simulation particles in each run")
    call CFG_add(cfg, "output%particles_dt", 1e-9_dp, &
         "Output time step for particles (s)")

  end subroutine pdsim_create_config

  !> Initialize pdsim module and the modules it depends on
  subroutine pdsim_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n, i

    ! Mesh related parameters
    character(len=200)              :: mesh_file
    real(dp)                        :: r_scale_factor, E_scale_factor
    character(len=120)              :: material_name
    character(len=120), allocatable :: field_component_names(:)
    integer                         :: n_field_comp

    call CFG_get(cfg, "output%name", pdsim_output_name)
    call CFG_get(cfg, "output%verbosity", pdsim_verbosity)
    call CFG_get(cfg, "output%level", pdsim_output_level)

    if (pdsim_output_level > 0) then
       call check_path_writable(trim(pdsim_output_name))
       call CFG_write(cfg, trim(pdsim_output_name) // "_config.cfg")
    end if

    call CFG_get(cfg, "input%mesh", mesh_file)
    call CFG_get(cfg, "input%coordinate_scale_factor", r_scale_factor)

    call iu_read_grid(trim(mesh_file), pdsim_ug, r_scale_factor)
    pdsim_ndim = iu_ndim_cell_type(pdsim_ug%cell_type)

    call CFG_get(cfg, "input%material_name", material_name)
    call CFG_get(cfg, "input%lookup_table_size", pdsim_table_size)

    call store_material_data(pdsim_ug, trim(material_name), &
         pdsim_icdata_material)

    ! Determine the number of materials, which should be numbered 0 (gas) to
    ! N, where N is the number of non-gas materials.
    pdsim_num_materials = 1 + &
         maxval(pdsim_ug%icell_data(:, pdsim_icdata_material))

    call CFG_get_size(cfg, "input%field_component_names", n_field_comp)

    if (pdsim_ndim /= n_field_comp) then
       write(error_unit, *) "Number of electric field components: ", &
            n_field_comp
       write(error_unit, *) "Problems dimension: ", pdsim_ndim
       error stop "Number of E-components must match dimension"
    end if

    allocate(field_component_names(n_field_comp))
    call CFG_get(cfg, "input%field_component_names", field_component_names)

    do n = 1, n_field_comp
       call iu_get_point_data_index(pdsim_ug, trim(field_component_names(n)), &
            pdsim_pdata_field(n))
       if (pdsim_pdata_field(n) == -1) then
          write(error_unit, *) trim(field_component_names(n)) // " not found"
          write(error_unit, *) "Available are:"
          do i = 1, pdsim_ug%n_point_data
             write(error_unit, *) i, trim(pdsim_ug%point_data_names(i))
          end do
          error stop "invalid input%field_component_names"
       end if
    end do

    call CFG_get(cfg, "input%field_scale_factor", E_scale_factor)
    if (abs(E_scale_factor - 1.0_dp) > 0) then
       pdsim_ug%point_data(:, pdsim_pdata_field(1:n_field_comp)) = E_scale_factor * &
            pdsim_ug%point_data(:, pdsim_pdata_field(1:n_field_comp))
    end if

    call CFG_get(cfg, "input%axisymmetric", pdsim_axisymmetric)

    if (pdsim_ndim == 2) then
       if (pdsim_axisymmetric) then
          pdsim_coord_system = pdsim_coord_axi
       else
          pdsim_coord_system = pdsim_coord_2d
       end if
    else
       pdsim_coord_system = pdsim_coord_3d
    end if

    call CFG_get_size(cfg, "input%ion_gamma", n)

    if (n /= 0 .and. n /= pdsim_num_materials) then
       print *, "input%ion_gamma should be specified for each"
       print *, "material (domain boundary, material 1, ..., material N)"
       print *, "Got size ", n, " while expecting size ", pdsim_num_materials
       error stop "Invalid size for input%ion_gamma"
    end if

    allocate(pdsim_ion_gamma(n))
    call CFG_get(cfg, "input%ion_gamma", pdsim_ion_gamma)
    pdsim_ion_see_enabled = (maxval(pdsim_ion_gamma) > 0.0_dp)

    call CFG_get_size(cfg, "input%photoemission_gamma", n)

    if (n /= 0 .and. n /= pdsim_num_materials) then
       print *, "input%photoemission_gamma should be specified for each"
       print *, "material (domain boundary, material 1, ..., material N)"
       print *, "Got size ", n, " while expecting size ", pdsim_num_materials
       error stop "Invalid size for input%photoemission_gamma"
    end if

    allocate(pdsim_photoemission_gamma(n))
    call CFG_get(cfg, "input%photoemission_gamma", pdsim_photoemission_gamma)
    pdsim_photoemission_enabled = (maxval(pdsim_photoemission_gamma) > 0.0_dp)

    if (maxval(pdsim_photoemission_gamma) > 1) &
         error stop "input%photoemission_gamma should be <= 1.0"

    call read_transport_data(cfg)

  end subroutine pdsim_initialize

  !> Read input data and store it in a lookup table
  subroutine read_transport_data(cfg)
    type(CFG_t), intent(inout) :: cfg

    character(len=200)             :: td_file
    character(len=20)              :: interp_method, xspacing_method
    integer                        :: input_interpolation, n, n_third
    integer                        :: xspacing
    character(len=20)              :: first_species
    character(len=20), allocatable :: third_bodies(:)
    real(dp), allocatable          :: xx(:), yy(:), efficiencies(:), v_drift(:)
    real(dp)                       :: max_field, factor
    real(dp)                       :: temperature, pressure
    character(len=20), allocatable :: gas_components(:)
    real(dp), allocatable          :: gas_fracs(:)

    call CFG_get(cfg, "input%transport_data_file", td_file)
    call CFG_get(cfg, "input%input_interpolation", interp_method)
    call CFG_get(cfg, "input%interpolation_xspacing", xspacing_method)

    select case (interp_method)
    case ("linear")
       input_interpolation = table_interp_linear
    case ("cubic_spline")
       input_interpolation = table_interp_cubic_spline
    case default
       error stop "Invalid input_interpolation method"
    end select

    select case (xspacing_method)
    case ("linear")
       xspacing = LT_xspacing_linear
    case ("quadratic")
       xspacing = LT_xspacing_quadratic
    case ("cubic")
       xspacing = LT_xspacing_cubic
    case default
       error stop "Invalid input%interpolation_xspacing"
    end select

    ! Set the gas properties
    call CFG_get(cfg, "gas%temperature", temperature)
    call CFG_get(cfg, "gas%pressure", pressure)

    call table_strings_numbers_from_file(td_file, "gas_composition", &
         gas_components, gas_fracs)

    if (pdsim_verbosity > 0) then
       print *, "Gas composition:"
       print *, "----------------------------"
       do n = 1, size(gas_components)
          write(*, "(A,A20,F8.3)") " ", gas_components(n), gas_fracs(n)
       end do
       print *, "----------------------------"
    end if

    call GAS_initialize(gas_components, gas_fracs, pressure, temperature)

    if (pdsim_verbosity > 0) then
       write(*, "(A,E11.3)") " Gas number density (1/m3): ", GAS_number_dens
    end if

    ! Set ion mobility
    call CFG_get(cfg, "input%ion_mobility", pdsim_ion_mobility)

    ! Convert to non-reduced units
    pdsim_ion_mobility = pdsim_ion_mobility / GAS_number_dens

    ! Start reading transport data
    call table_from_file(td_file, "Mobility *N (1/m/V/s)", xx, yy)

    xx = xx * Townsend_to_SI * GAS_number_dens
    max_field = xx(size(xx))
    yy = yy / GAS_number_dens
    pdsim_tdtbl = LT_create(0.0_dp, max_field, pdsim_table_size, pdsim_ncols, &
         xspacing)
    call table_set_column(pdsim_tdtbl, pdsim_col_mu, xx, yy, &
         input_interpolation)

    call table_from_file(td_file, "Townsend ioniz. coef. alpha/N (m2)", &
         xx, yy)
    xx = xx * Townsend_to_SI * GAS_number_dens
    yy = yy * GAS_number_dens
    call table_set_column(pdsim_tdtbl, pdsim_col_alpha, xx, yy, &
         input_interpolation)

    call table_from_file(td_file, "Townsend attach. coef. eta/N (m2)", &
         xx, yy)
    xx = xx * Townsend_to_SI * GAS_number_dens
    yy = yy * GAS_number_dens
    call table_set_column(pdsim_tdtbl, pdsim_col_eta, xx, yy, &
         input_interpolation)

    ! Check for three-body attachment
    call CFG_get(cfg, "input%three_body_first_species", first_species)
    call CFG_get_size(cfg, "input%three_body_second_species", n_third)

    if (n_third > 0) then
       allocate(third_bodies(n_third))
       allocate(efficiencies(n_third))
       call CFG_get(cfg, "input%three_body_second_species", third_bodies)
       call CFG_get(cfg, "input%three_body_efficiencies", efficiencies)

       call table_from_file(td_file, "Three-body attachment rate (m6/s)", &
            xx, yy)

       xx = xx * Townsend_to_SI * GAS_number_dens
       v_drift = xx * LT_get_col(pdsim_tdtbl, pdsim_col_mu, xx)

       ! Convert rate to Townsend attachment coefficient
       yy = yy / v_drift

       ! Multiply with number densities and efficiencies of third bodies
       factor = 0.0_dp
       do n = 1, n_third
          factor = factor + GAS_number_dens * &
               GAS_get_fraction(third_bodies(n)) * efficiencies(n)
       end do

       ! Multiply with first species
       factor = factor * GAS_number_dens * GAS_get_fraction(first_species)
       yy = yy * factor

       call table_set_column(pdsim_tdtbl, pdsim_col_eta, xx, yy, &
            input_interpolation, add=.true.)
    end if

    if (pdsim_output_level > 0) then
       call write_transport_data_summary(&
            trim(pdsim_output_name) // "_transport_data.txt")
    end if

  end subroutine read_transport_data

  !> Write the transport data as used by the code to a file
  subroutine write_transport_data_summary(fname)
    character(len=*), intent(in) :: fname
    integer                      :: my_unit, n

    open(newunit=my_unit, file=trim(fname), action="write")
    write(my_unit, "(A)") "E[V/m] alpha[1/m] eta[1/m] mu[m^2/(Vs)]"
    do n = 1, pdsim_tdtbl%n_points
       write(my_unit, *) pdsim_tdtbl%x(n), &
            pdsim_tdtbl%cols_rows(pdsim_col_alpha, n), &
            pdsim_tdtbl%cols_rows(pdsim_col_eta, n), &
            pdsim_tdtbl%cols_rows(pdsim_col_mu, n)
    end do
    write(my_unit, *) ""
    close(my_unit)

    if (pdsim_verbosity > 0) print *, "Wrote ", trim(fname)

    do n = 1, pdsim_tdtbl%n_points
       if (pdsim_tdtbl%cols_rows(pdsim_col_alpha, n) > &
            pdsim_tdtbl%cols_rows(pdsim_col_eta, n)) exit
    end do

    if (n <= pdsim_tdtbl%n_points) then
       if (pdsim_verbosity > 0) then
          write(*, "(A,2E11.3)") " Critical field (V/m): ", pdsim_tdtbl%x(n)
       end if
    else
       error stop "Critical not reached in input data"
    end if
  end subroutine write_transport_data_summary

  subroutine check_path_writable(pathname)
    character(len=*), intent(in) :: pathname
    integer                      :: my_unit, iostate

    open(newunit=my_unit, file=trim(pathname)//"_DUMMY", iostat=iostate)
    if (iostate /= 0) then
       print *, "Output file: " // trim(pathname)
       error stop "Directory not writable (does it exist?)"
    else
       close(my_unit, status='delete')
    end if
  end subroutine check_path_writable

  !> Convert particle coordinates to (r, z)
  pure function x_to_rz(x) result(rz)
    real(dp), intent(in) :: x(3)
    real(dp)             :: rz(2)

    rz(1) = sqrt(x(1)**2 + x(3)**2) ! radius
    rz(2) = x(2)                    ! z
  end function x_to_rz

  !> Convert 3D coordinate xyz to used actual coordinate system
  pure function pdsim_convert_r(xyz) result(r)
    real(dp), intent(in) :: xyz(3)
    real(dp)             :: r(3)

    select case (pdsim_coord_system)
    case (pdsim_coord_2d)
       r(1:2) = xyz(1:2)
       r(3) = 0.0_dp
    case (pdsim_coord_axi)
       r(1:2) = x_to_rz(xyz)
       r(3) = 0.0_dp
    case (pdsim_coord_3d)
       r = xyz
    end select
  end function pdsim_convert_r

  !> If the 'material' variable is not "NONE", find the cell data or point
  !> data corresponding to it, and store it as integer cell data.
  subroutine store_material_data(ug, name, i_icdata)
    type(iu_grid_t), intent(inout)  :: ug
    !> Name of material variable or "NONE" if the whole domain is gas
    character(len=*), intent(in)    :: name
    !> Index of integer cell data variable corresponding to material
    integer, intent(out)            :: i_icdata
    real(dp)                        :: fac
    integer                         :: n, i_cdata, i_pdata

    call iu_get_icell_data_index(ug, trim(name), i_icdata)
    call iu_get_cell_data_index(ug, trim(name), i_cdata)
    call iu_get_point_data_index(ug, trim(name), i_pdata)

    if (name == "NONE") then
       ! Add all-gas material
       call iu_add_icell_data(ug, "material", i_icdata)
       ug%icell_data(:, i_icdata) = pdsim_gas_material_value
    else if (i_icdata == -1 .and. i_cdata > 0) then
       ! Convert cell data to integer
       call iu_add_icell_data(ug, "material", i_icdata)
       ug%icell_data(:, i_icdata) = nint(ug%cell_data(:, i_cdata))
    else if (i_icdata == -1 .and. i_pdata > 0) then
       ! Convert point data to cell data
       call iu_add_icell_data(ug, "material", i_icdata)

       fac = 1.0_dp/ug%n_points_per_cell
       do n = 1, ug%n_cells
          ug%icell_data(n, i_icdata) = nint(&
               sum(ug%point_data(ug%cells(:, n), i_pdata)) * fac)
       end do
    else if (i_icdata == -1) then
       write(error_unit, *) trim(name) // " not found"
       error stop "invalid input%material_name"
    end if

  end subroutine store_material_data

  !> Compute volume average of a variable defined at points (vertices), only
  !> considering the gas phase
  subroutine pdsim_pointdata_average(ug, iv, axisymmetric, power, avg, vol)
    type(iu_grid_t), intent(in) :: ug
    integer, intent(in)         :: iv    !< Index of point data variable
    logical, intent(in)         :: axisymmetric
    integer, intent(in)         :: power !< Use dV**power (1 for normal average)
    real(dp), intent(out)       :: avg   !< Average of the variable
    real(dp), intent(out)       :: vol   !< Total volume
    integer                     :: n, k, i_point
    real(dp)                    :: dV, total_sum, w
    real(dp)                    :: center(3)
    real(dp), parameter         :: pi = acos(-1.0_dp)
    logical, allocatable        :: mask(:)

    if (axisymmetric .and. iu_ndim_cell_type(ug%cell_type) /= 2) &
         error stop "Axisymmetric should only be used in 2D"

    allocate(mask(ug%n_cells))
    mask(:) = (ug%icell_data(:, pdsim_icdata_material) == &
         pdsim_gas_material_value)

    ! Weight per point. An improvement could be to have a variable weight in
    ! axisymmetric cases.
    w = 1.0_dp / ug%n_points_per_cell

    vol = 0.0_dp
    total_sum = 0.0_dp

    do n = 1, ug%n_cells
       if (mask(n)) then
          if (axisymmetric) then
             ! Multiply with 2 * pi * r, where r is at the cell center
             center = iu_get_cell_center(ug, n)
             dV = ug%cell_volume(n) * 2 * pi * center(1)
          else
             dV = ug%cell_volume(n)
          end if

          vol = vol + dV
          do k = 1, ug%n_points_per_cell
             i_point = ug%cells(k, n)
             total_sum = total_sum + dV**power * w * &
                  ug%point_data(i_point, iv)
          end do
       end if
    end do

    avg = total_sum / vol**power

  end subroutine pdsim_pointdata_average

end module m_pdsim
