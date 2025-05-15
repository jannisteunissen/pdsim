module m_pdsim
  use iso_fortran_env, only: error_unit
  use m_interp_unstructured
  use m_particle_core
  use m_config
  use m_photoi
  use m_units_constants
  use m_cross_sec
  use m_gas

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  real(dp), parameter :: undefined_real = -1e100_dp
  character(len=*), parameter :: undefined_str = "NOT_SPECIFIED"

  integer, parameter :: pdsim_coord_2d = 1
  integer, parameter :: pdsim_coord_3d = 2
  integer, parameter :: pdsim_coord_axi = 3

  real(dp), parameter :: material_threshold = 1e-8_dp

  type pdsim_t
     integer  :: n_initial_positions
     real(dp) :: initial_rmin(3)
     real(dp) :: initial_rmax(3)

     real(dp) :: max_dt
     real(dp) :: num_electrons_inception

     character(len=200) :: output_name
     logical            :: output_log
     logical            :: output_particles
     real(dp)           :: output_particles_dt

     type(PC_t) :: pc
  end type pdsim_t

  type(iu_grid_t) :: pdsim_ug
  integer         :: i_pdata_field(3)
  integer         :: i_cdata_material
  integer         :: pdsim_gas_material_value
  integer         :: pdsim_coord_system
  integer         :: pdsim_ndim
  logical         :: pdsim_axisymmetric

  ! Public types
  public :: dp
  public :: pdsim_t

  ! Public data
  public :: pdsim_ug

  ! Public methods
  public :: pdsim_create_config
  public :: pdsim_initialize
  public :: pdsim_get_gas_cells
  public :: pdsim_write_particles
  public :: pdsim_handle_events

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
         "Value of the material variable in the gas phase")
    call CFG_add(cfg, "input%lookup_table_size", 1000, &
         "Size to use for cross section lookup table")

    call CFG_add(cfg, "electrons%max_eV", 500.0_dp, &
         "Max. electron energy (eV)")

    call CFG_add(cfg, "gas%temperature", 300.0_dp, "Gas temperature (K)")
    call CFG_add(cfg, "gas%pressure", 1.0_dp, "Gas pressure (bar)")
    call CFG_add(cfg, "gas%components", ["N2", "O2"], "Gas components", .true.)
    call CFG_add(cfg, "gas%fractions", [0.8_dp, 0.2_dp], &
         "Partial pressure of the gases (as if they were ideal gases)", .true.)
    call CFG_add(cfg, "gas%cross_section_file", undefined_str, &
         "File containing electron-neutral cross sections", required=.true.)

    call CFG_add(cfg, "output%name", "output/pdsim", &
         "The base name for output files")
    call CFG_add(cfg, "output%log", .true., &
         "Write output to a log file")
    call CFG_add(cfg, "output%particles", .true., &
         "Output the simulation particles in each run")
    call CFG_add(cfg, "output%particles_dt", 1e-9_dp, &
         "Output time step for particles (s)")

    call CFG_add(cfg, "simulation%max_dt", 1.0e-12_dp, &
         "Maximal time step (s)")
    call CFG_add(cfg, "simulation%n_initial_positions", 1, &
         "Number of initial electron positions to consider")
    call CFG_add(cfg, "simulation%num_electrons_inception", 1e5_dp, &
         "Assume inception takes place when there are this many electrons")
    call CFG_add(cfg, "simulation%initial_rmin", &
         [undefined_real, undefined_real, undefined_real], &
         "Min. coordinate for generating initial electrons")
    call CFG_add(cfg, "simulation%initial_rmax", &
         [undefined_real, undefined_real, undefined_real], &
         "Max. coordinate for generating initial electrons")

    call photoi_create_cfg(cfg)

  end subroutine pdsim_create_config

  !> Initialize pdsim module and the modules it depends on
  subroutine pdsim_initialize(cfg, pd)
    type(CFG_t), intent(inout)   :: cfg
    type(pdsim_t), intent(inout) :: pd

    integer :: n, i

    ! Mesh related parameters
    character(len=200) :: mesh_file
    real(dp)           :: r_scale_factor, E_scale_factor
    character(len=80)  :: material_name
    character(len=80), allocatable :: field_component_names(:)
    integer            :: n_field_comp

    ! Gas related parameters
    integer                        :: n_gas_comp, n_gas_frac
    character(len=20), allocatable :: gas_names(:)
    real(dp), allocatable          :: gas_fracs(:)
    real(dp)                       :: temperature, pressure
    character(len=200)             :: cs_file
    type(CS_t), allocatable        :: cross_secs(:)

    ! Particle model related parameters
    real(dp) :: max_eV
    integer :: lookup_table_size

    call CFG_get(cfg, "input%mesh", mesh_file)
    call CFG_get(cfg, "input%coordinate_scale_factor", r_scale_factor)

    call iu_read_grid(trim(mesh_file), pdsim_ug, r_scale_factor)
    pdsim_ndim = iu_ndim_cell_type(pdsim_ug%cell_type)

    call CFG_get(cfg, "input%gas_material_value", pdsim_gas_material_value)
    call CFG_get(cfg, "input%material_name", material_name)

    call store_material_as_cell_data(pdsim_ug, trim(material_name), &
            i_cdata_material)

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
            i_pdata_field(n))
       if (i_pdata_field(n) == -1) then
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
       pdsim_ug%point_data(:, i_pdata_field(1:n_field_comp)) = E_scale_factor * &
             pdsim_ug%point_data(:, i_pdata_field(1:n_field_comp))
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
    call CFG_get(cfg, "gas%cross_section_file", cs_file)

    call GAS_initialize(gas_names, gas_fracs, pressure, temperature)

    ! Load cross sections
    do n = 1, GAS_num_gases
       call CS_add_from_file(trim(cs_file), trim(gas_names(n)), &
            gas_fracs(n) * GAS_number_dens, 0.0_dp, cross_secs)
    end do

    call CFG_get(cfg, "simulation%initial_rmin", pd%initial_rmin)
    call CFG_get(cfg, "simulation%initial_rmax", pd%initial_rmax)

    if (all(pd%initial_rmin <= undefined_real)) &
         pd%initial_rmin = pdsim_ug%rmin
    if (all(pd%initial_rmax <= undefined_real)) &
         pd%initial_rmax = pdsim_ug%rmax

    call CFG_get(cfg, "electrons%max_eV", max_eV)

    call CFG_get(cfg, "input%lookup_table_size", lookup_table_size)
    call CFG_get(cfg, "simulation%num_electrons_inception", &
         pd%num_electrons_inception)
    call CFG_get(cfg, "simulation%max_dt", pd%max_dt)
    call CFG_get(cfg, "simulation%n_initial_positions", pd%n_initial_positions)

    call CFG_get(cfg, "output%name", pd%output_name)
    call check_path_writable(trim(pd%output_name))

    call CFG_get(cfg, "output%log", pd%output_log)
    call CFG_get(cfg, "output%particles", pd%output_particles)
    call CFG_get(cfg, "output%particles_dt", pd%output_particles_dt)

    if (pdsim_axisymmetric) then
       pd%pc%particle_mover => PC_verlet_cyl_advance
       pd%pc%after_mover    => PC_verlet_cyl_correct_accel
    else
       pd%pc%particle_mover => PC_verlet_advance
       pd%pc%after_mover    => PC_verlet_correct_accel
    end if

    pd%pc%accel_function => accel_function
    pd%pc%outside_check => outside_check

    call pd%pc%initialize(UC_elec_mass, 100*int(pd%num_electrons_inception))
    call pd%pc%use_cross_secs(max_eV, lookup_table_size, cross_secs)

    where (pd%pc%colls(:)%type == CS_ionize_t .or. &
         pd%pc%colls(:)%type == CS_attach_t)
       pd%pc%coll_is_event(:) = .true.
    end where

    call photoi_initialize(cfg)

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

  !> Write the particle coordinates to a file that can be loaded in Visit
  subroutine pdsim_write_particles(pd, i_step, i_run)
    type(pdsim_t), intent(in) :: pd
    integer, intent(in)       :: i_step
    integer, intent(in)       :: i_run
    character(len=200)        :: fname
    integer                   :: my_unit, n

    write(fname, "(A,I0.5,A,I0.5,A)") trim(pd%output_name) // '_run', &
         i_run, "_", i_step, ".3D"

    open(newunit=my_unit, file=trim(fname), action='write')
    write(my_unit, *) "X Y Z eV"
    do n = 1, pd%pc%n_part
       write(my_unit, *) pd%pc%particles(n)%x, &
            PC_v_to_en(pd%pc%particles(n)%v, pd%pc%mass) / UC_elec_volt
    end do
    close(my_unit)

    print *, "Wrote ", trim(fname)
  end subroutine pdsim_write_particles

  !> Construct index array of approximately n_max cells that are in the gas
  subroutine pdsim_get_gas_cells(ug, n_max, rmin, rmax, i_cell_gas)
    type(iu_grid_t)                     :: ug
    integer, intent(in)                 :: n_max
    real(dp), intent(in)                :: rmin(3)
    real(dp), intent(in)                :: rmax(3)
    integer, allocatable, intent(inout) :: i_cell_gas(:)
    logical, allocatable                :: mask(:)
    real(dp), allocatable               :: unif_01(:)
    integer, allocatable                :: indices(:)
    integer                             :: n, i, n_gas_cells
    real(dp)                            :: center(3)

    if (i_cdata_material > 0) then
       mask = abs(ug%cell_data(:, i_cdata_material) - &
            pdsim_gas_material_value) < material_threshold
    else
       allocate(mask(ug%n_cells))
       mask(:) = .true.
    end if

    do n = 1, ug%n_cells
       center = iu_get_cell_center(ug, n)
       if (any(center < rmin) .or. any(center > rmax)) then
          mask(n) = .false.
       end if
    end do

    n_gas_cells = count(mask)
    allocate(i_cell_gas(n_gas_cells))

    i = 0
    do n = 1, ug%n_cells
       if (mask(n)) then
          i = i + 1
          i_cell_gas(i) = n
       end if
    end do

    if (n_gas_cells > n_max) then
       ! Randomly sub-sample
       allocate(unif_01(n_max))
       call random_number(unif_01)

       ! unif_01 are in the range [0, 1)
       indices = 1 + floor(unif_01 * n_gas_cells)
       i_cell_gas = i_cell_gas(indices)
    else if (n_gas_cells == 0) then
       error stop "No valid initial cells"
    end if

  end subroutine pdsim_get_gas_cells

  !> Check whether a particle is outside the gas
  integer function outside_check(my_part)
    type(PC_part_t), intent(inout) :: my_part
    integer                        :: i_cell
    real(dp)                       :: material, x(3)

    outside_check = 0

    if (i_cdata_material > 0) then
       select case (pdsim_coord_system)
       case (pdsim_coord_2d)
          x(1:2) = my_part%x(1:2)
          x(3) = 0.0_dp
       case (pdsim_coord_axi)
          x(1:2) = x_to_rz(my_part%x)
          x(3) = 0.0_dp
       case (pdsim_coord_3d)
          x = my_part%x
       end select

       i_cell = 0
       material = -1e100_dp
       call iu_get_cell_scalar_at(pdsim_ug, x, i_cdata_material, &
            i_cell, material)

       if (abs(material - pdsim_gas_material_value) > 1e-8_dp) then
          outside_check = 1
       end if
    end if

  end function outside_check

  !> Get acceleration of particle by interpolating the electric field
  function accel_function(my_part) result(a)
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: a(3), x(3)

    select case (pdsim_coord_system)
    case (pdsim_coord_2d)
       x(1:2) = my_part%x(1:2)
       x(3) = 0.0_dp
       a(3) = 0.0_dp
    case (pdsim_coord_axi)
       x(1:2) = x_to_rz(my_part%x)
       x(3) = 0.0_dp
       a(3) = 0.0_dp
    case (pdsim_coord_3d)
       x = my_part%x
    end select

    call iu_interpolate_at(pdsim_ug, x, pdsim_ndim, &
         i_pdata_field(1:pdsim_ndim), a(1:pdsim_ndim), my_part%id)
    a = a * UC_elec_q_over_m

    if (my_part%id < 1) error stop "accel_function: interpolation error"
  end function accel_function

  !> Convert particle coordinates to (r, z)
  pure function x_to_rz(x) result(rz)
    real(dp), intent(in) :: x(3)
    real(dp)             :: rz(2)

    rz(1) = sqrt(x(1)**2 + x(3)**2) ! radius
    rz(2) = x(2)                    ! z
  end function x_to_rz

  !> After updating the particles, events such as ionization are stored. Here
  !> we compute photoionization based on these events.
  subroutine pdsim_handle_events(pc, n_ionizations)
    type(PC_t), intent(inout) :: pc
    integer, intent(inout)    :: n_ionizations

    integer                     :: n, n_photons
    real(dp), allocatable, save :: coords(:, :)
    real(dp), allocatable, save :: weights(:)
    type(PC_part_t)             :: my_part
    real(dp), parameter         :: array_incr_fac = 2.0_dp

    n = max(pc%n_events, pc%n_part)

    if (.not. allocated(weights)) then
       n = nint(n * array_incr_fac)
       allocate(coords(3, n))
       allocate(weights(n))
    else if (size(weights) < n) then
       n = nint(n * array_incr_fac)
       deallocate(coords)
       deallocate(weights)
       allocate(coords(3, n))
       allocate(weights(n))
    end if

    n_ionizations = n_ionizations + &
         count(pc%event_list(1:pc%n_events)%ctype == CS_ionize_t)

    if (photoi_enabled) then
       call photoi_from_events(pc%n_events, pc%event_list, pc%rng, &
            coords, weights, n_photons)

       do n = 1, n_photons
          if (pdsim_axisymmetric) then
             my_part%x(1) = norm2(coords(1:2, n))
             my_part%x(2) = coords(3, n)
             my_part%x(3) = 0.0_dp
          else
             my_part%x(1:pdsim_ndim) = coords(1:pdsim_ndim, n)
             my_part%x(pdsim_ndim+1:) = 0.0_dp
          end if

          if (outside_check(my_part) == 0) then
             my_part%v(:)   = 0.0_dp
             my_part%w      = 1.0_dp
             my_part%t_left = 0.0_dp
             my_part%a(:)   = 0.0_dp ! Will be set later
             call pc%add_part(my_part)
          end if
       end do
    end if

    ! do n = 1, pc%n_events
    !    if (pc%event_list(n)%ctype == CS_attach_t) then
    !       ! Create negative ion
    !       my_part = pc%event_list(n)%part
    !       my_part%v = ion_velocity(my_part)
    !       my_part%a = 0.0_dp
    !       call pc_ion%add_part(my_part)
    !    end if
    ! end do

    pc%n_events = 0

  end subroutine pdsim_handle_events

  !> If the 'material' variable is not "NONE", find the cell data or point
  !> data corresponding to it. In case of point data, convert to cell data.
  subroutine store_material_as_cell_data(ug, name, i_material)
    type(iu_grid_t), intent(inout)  :: ug
    !> Name of material variable or "NONE" if the whole domain is gas
    character(len=*), intent(in)    :: name
    !> Index of cell data variable corresponding to material
    integer, intent(out)            :: i_material
    real(dp), allocatable           :: old_data(:, :)
    character(len=128), allocatable :: old_names(:)
    real(dp)                        :: fac
    integer                         :: n

    if (name == "NONE") then
       i_material = -1
       return
    end if

    call iu_get_cell_data_index(ug, trim(name), i_material)
    if (i_material >= 1) return

    call iu_get_point_data_index(ug, trim(name), i_material)

    if (i_material == -1) then
       write(error_unit, *) trim(name) // " not found"
       error stop "invalid input%material_name"
    end if

    ! Convert point data to cell data
    call move_alloc(ug%cell_data, old_data)
    call move_alloc(ug%cell_data_names, old_names)
    allocate(ug%cell_data(ug%n_cells, ug%n_cell_data+1))
    allocate(ug%cell_data_names(ug%n_cell_data+1))

    ug%cell_data(:, 1:ug%n_cell_data) = old_data
    ug%cell_data_names(1:ug%n_cell_data) = old_names

    ug%n_cell_data = ug%n_cell_data + 1
    ug%cell_data_names(ug%n_cell_data) = name

    fac = 1.0_dp/ug%n_points_per_cell

    do n = 1, ug%n_cells
       ug%cell_data(n, ug%n_cell_data) = &
            sum(ug%point_data(ug%cells(:, n), i_material)) * fac
    end do

    i_material = ug%n_cell_data

  end subroutine store_material_as_cell_data

end module m_pdsim
