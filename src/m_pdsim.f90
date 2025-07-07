module m_pdsim
  use iso_fortran_env, only: error_unit, int64
  use m_interp_unstructured
  use m_particle_core
  use m_config
  use m_cross_sec
  use m_gas
  use m_lookup_table

  implicit none
  private

  integer, parameter, public :: dp = kind(0.0d0)
  real(dp), parameter, public :: undefined_real = -1e100_dp
  character(len=*), parameter, public :: undefined_str = "NOT_SPECIFIED"

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

  integer, public :: i_k_integral, i_alpha, i_eta
  integer, public :: i_w, i_p0
  integer, public :: i_x1, i_x2, i_x3
  integer, public :: i_ion_x1, i_ion_x2, i_ion_x3
  integer, public :: i_avalanche_time
  integer, public :: i_ion_time
  integer, public :: i_inception_time
  integer, public :: i_inception_prob

  character(len=200), public, protected :: pdsim_output_name

  !> Size for lookup tables used for interpolation
  integer, public, protected :: pdsim_table_size

  !> Positive ion mobility (m^2/Vs)
  real(dp), public, protected :: pdsim_ion_mobility

  !> Secondary emission coefficient for ions reaching a boundary
  real(dp), public, protected :: pdsim_ion_gamma_boundary

  type(iu_grid_t), public    :: pdsim_ug
  integer, public, protected :: pdsim_pdata_field(3)
  integer, public, protected :: pdsim_cdata_material
  integer, public, protected :: pdsim_pdata_material
  integer, public, protected :: pdsim_gas_material_value
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
    call CFG_add(cfg, "input%gas_material_value", 0, &
         "Value of the 'gas' material, should be lower than other materials")
    call CFG_add(cfg, "input%lookup_table_size", 1000, &
         "Size to use for cross section lookup table")
    call CFG_add(cfg, "input%alpha_eta_file", undefined_str, &
         "File with field (V/m), alpha (1/m), eta (1/m), mu (m^2/Vs)")
    call CFG_add(cfg, "input%ion_mobility", 2e-4_dp, &
         "Positive ion mobility (m^2/Vs)")
    call CFG_add(cfg, "input%ion_gamma_boundary", 0.0_dp, &
         "Secondary emission coefficient for ions reaching a boundary")

    call CFG_add(cfg, "gas%temperature", 300.0_dp, "Gas temperature (K)")
    call CFG_add(cfg, "gas%pressure", 1.0_dp, "Gas pressure (bar)")
    call CFG_add(cfg, "gas%components", ["N2", "O2"], "Gas components", .true.)
    call CFG_add(cfg, "gas%fractions", [0.8_dp, 0.2_dp], &
         "Partial pressure of the gases (as if they were ideal gases)", .true.)

    call CFG_add(cfg, "output%name", "output/pdsim", &
         "The base name for output files")
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

    ! Gas related parameters
    integer                        :: n_gas_comp, n_gas_frac
    character(len=20), allocatable :: gas_names(:)
    real(dp), allocatable          :: gas_fracs(:)
    real(dp)                       :: temperature, pressure

    ! Transport data
    character(len=200)    :: alpha_eta_file
    integer               :: n_rows
    real(dp), allocatable :: field_alpha_eta(:, :)

    call CFG_get(cfg, "output%name", pdsim_output_name)
    call check_path_writable(trim(pdsim_output_name))

    call CFG_write(cfg, trim(pdsim_output_name) // "_config.cfg")

    call CFG_get(cfg, "input%mesh", mesh_file)
    call CFG_get(cfg, "input%coordinate_scale_factor", r_scale_factor)

    call iu_read_grid(trim(mesh_file), pdsim_ug, r_scale_factor)
    pdsim_ndim = iu_ndim_cell_type(pdsim_ug%cell_type)

    call CFG_get(cfg, "input%gas_material_value", pdsim_gas_material_value)
    call CFG_get(cfg, "input%material_name", material_name)
    call CFG_get(cfg, "input%lookup_table_size", pdsim_table_size)

    call store_material_data(pdsim_ug, trim(material_name), &
         pdsim_cdata_material, pdsim_pdata_material)

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

    call CFG_get(cfg, "input%ion_mobility", pdsim_ion_mobility)
    call CFG_get(cfg, "input%ion_gamma_boundary", pdsim_ion_gamma_boundary)

    call CFG_get_size(cfg, "gas%components", n_gas_comp)
    call CFG_get_size(cfg, "gas%fractions", n_gas_frac)
    if (n_gas_comp /= n_gas_frac) &
         error stop "gas%components and gas%fractions have unequal size"
    allocate(gas_names(n_gas_comp))
    allocate(gas_fracs(n_gas_comp))

    call CFG_get(cfg, "gas%temperature", temperature)
    call CFG_get(cfg, "gas%pressure", pressure)
    call CFG_get(cfg, "gas%components", gas_names)
    call CFG_get(cfg, "gas%fractions", gas_fracs)

    call GAS_initialize(gas_names, gas_fracs, pressure, temperature)

    call CFG_get(cfg, "input%alpha_eta_file", alpha_eta_file)

    if (alpha_eta_file == undefined_str) &
         error stop "input%alpha_eta_file is not set"

    call read_table_from_txt(trim(alpha_eta_file), pdsim_ncols+1, 1000, &
         field_alpha_eta)
    n_rows = size(field_alpha_eta, 2)

    pdsim_tdtbl = LT_create(field_alpha_eta(1, 1), &
         field_alpha_eta(1, n_rows), pdsim_table_size, pdsim_ncols)

    call LT_set_col(pdsim_tdtbl, pdsim_col_alpha, field_alpha_eta(1, :), &
         field_alpha_eta(2, :))
    call LT_set_col(pdsim_tdtbl, pdsim_col_eta, field_alpha_eta(1, :), &
         field_alpha_eta(3, :))
    call LT_set_col(pdsim_tdtbl, pdsim_col_mu, field_alpha_eta(1, :), &
         field_alpha_eta(4, :))

  end subroutine pdsim_initialize

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
  !> data corresponding to it. Convert one from the other so both are present.
  subroutine store_material_data(ug, name, i_cdata, i_pdata)
    type(iu_grid_t), intent(inout)  :: ug
    !> Name of material variable or "NONE" if the whole domain is gas
    character(len=*), intent(in)    :: name
    !> Index of cell data variable corresponding to material
    integer, intent(out)            :: i_cdata
    !> Index of point data variable corresponding to material
    integer, intent(out)            :: i_pdata
    real(dp)                        :: fac
    integer                         :: n, k, i_point

    call iu_get_cell_data_index(ug, trim(name), i_cdata)
    call iu_get_point_data_index(ug, trim(name), i_pdata)

    if (name == "NONE") then
       ! Add all-gas material
       call iu_add_cell_data(ug, "material_cdata", i_cdata)
       call iu_add_point_data(ug, "material_pdata", i_pdata)

       ug%cell_data(:, i_cdata) = pdsim_gas_material_value
       ug%cell_data(:, i_pdata) = pdsim_gas_material_value

    else if (i_cdata > 0 .and. i_pdata == -1) then
       ! Add point data
       call iu_add_point_data(ug, "material_pdata", i_pdata)


       ! Set points to the lowest material value of the cells they are part of
       ug%point_data(:, i_pdata) = 1000*1000

       do n = 1, ug%n_cells
          do k = 1, ug%n_points_per_cell
             i_point = ug%cells(k, n)
             ug%point_data(i_point, i_pdata) = min(ug%cell_data(n, i_cdata), &
                  ug%point_data(i_point, i_pdata))
          end do
       end do

    else if (i_cdata == -1 .and. i_pdata > 0) then
       ! Add cell data
       call iu_add_cell_data(ug, "material_cdata", i_cdata)

       fac = 1.0_dp/ug%n_points_per_cell
       do n = 1, ug%n_cells
          ug%cell_data(n, ug%n_cell_data) = &
               sum(ug%point_data(ug%cells(:, n), i_cdata)) * fac
       end do

    else if (i_cdata == -1 .and. i_pdata == -1) then
       write(error_unit, *) trim(name) // " not found"
       error stop "invalid input%material_name"
    end if

  end subroutine store_material_data

  !> Routine to read in tabulated data from a text file
  subroutine read_table_from_txt(file_name, n_columns, max_rows, array)
    character(len=*), intent(in)       :: file_name
    integer, intent(in)                :: n_columns
    integer, intent(in)                :: max_rows

    real(dp), allocatable, intent(out) :: array(:, :)
    real(dp), allocatable              :: tmp_array(:, :)
    integer                            :: my_unit, n, io

    allocate(tmp_array(n_columns, max_rows))
    open(newunit=my_unit, file=trim(file_name), action = "read")

    do n = 1, max_rows
       read(my_unit, *, iostat=io) tmp_array(:, n)
       if (io /= 0) exit
    end do

    close(my_unit)

    allocate(array(n_columns, n-1))
    array(:, 1:n-1) = tmp_array(:, 1:n-1)
  end subroutine read_table_from_txt

  !> Compute volume average of a variable defined at points (vertices), only
  !> considering the gas phase
  subroutine pdsim_pointdata_average(ug, iv, axisymmetric, avg, vol)
    type(iu_grid_t), intent(in) :: ug
    integer, intent(in)         :: iv !< Index of point data variable
    logical, intent(in)         :: axisymmetric
    real(dp), intent(out)       :: avg !< Average of the variable
    real(dp), intent(out)       :: vol !< Total volume
    integer                     :: n, k, i_point
    real(dp)                    :: dV, total_sum, w
    real(dp)                    :: center(3)
    real(dp), parameter         :: pi = acos(-1.0_dp)
    logical, allocatable        :: mask(:)

    if (axisymmetric .and. iu_ndim_cell_type(ug%cell_type) /= 2) &
         error stop "Axisymmetric should only be used in 2D"

    allocate(mask(ug%n_cells))
    mask(:) = (nint(ug%cell_data(:, pdsim_cdata_material)) == &
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
             total_sum = total_sum + dV * w * ug%point_data(i_point, iv)
          end do
       end if
    end do

    avg = total_sum / vol

  end subroutine pdsim_pointdata_average

end module m_pdsim
