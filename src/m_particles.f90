!> Module that contains methods and data for particle-based simulations
module m_particles
  use m_config
  use m_interp_unstructured
  use m_pdsim
  use m_cross_sec
  use m_particle_core
  use m_units_constants
  use m_photoi
  use m_gas

  implicit none
  private

  logical, public, protected :: particles_enabled

  real(dp) :: max_dt
  real(dp) :: num_electrons_inception

  logical  :: particles_output
  real(dp) :: particles_output_dt

  type(PC_t) :: pc

  public :: particles_create_config
  public :: particles_initialize
  public :: particles_simulate

contains

  subroutine particles_create_config(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "particles%max_dt", 1.0e-12_dp, &
         "Maximal time step (s) for particles")
    call CFG_add(cfg, "particles%initial_position", [0.0_dp, 0.0_dp, 0.0_dp], &
         "Initial electron position (m); in 2D, third coordinate is zero)")
    call CFG_add(cfg, "particles%initial_energy", 0.0_dp, &
         "Initial energy of the first electron (eV)")
    call CFG_add(cfg, "particles%n_runs", 1, &
         "Number of runs")
    call CFG_add(cfg, "particles%num_electrons_inception", 1e5_dp, &
         "Assume inception takes place when there are this many electrons")
    call CFG_add(cfg, "particles%max_energy", 500.0_dp, &
         "Max. electron energy (eV)")
    call CFG_add(cfg, "particles%cross_section_file", undefined_str, &
         "File containing electron-neutral cross sections", required=.true.)

  end subroutine particles_create_config

  subroutine particles_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg

    integer                 :: n
    character(len=200)      :: cs_file
    type(CS_t), allocatable :: cross_secs(:)
    real(dp)                :: max_energy

    call CFG_get(cfg, "particles%cross_section_file", cs_file)

    ! Load cross sections
    do n = 1, GAS_num_gases
       call CS_add_from_file(trim(cs_file), trim(GAS_comp_names(n)), &
            GAS_comp_fracs(n) * GAS_number_dens, 0.0_dp, cross_secs)
    end do

    call CFG_get(cfg, "particles%max_energy", max_energy)

    call CFG_get(cfg, "particles%num_electrons_inception", &
         num_electrons_inception)
    call CFG_get(cfg, "particles%max_dt", max_dt)

    call CFG_get(cfg, "output%particles", particles_output)
    call CFG_get(cfg, "output%particles_dt", particles_output_dt)

    if (pdsim_axisymmetric) then
       pc%particle_mover => PC_verlet_cyl_advance
       pc%after_mover    => PC_verlet_cyl_correct_accel
    else
       pc%particle_mover => PC_verlet_advance
       pc%after_mover    => PC_verlet_correct_accel
    end if

    pc%accel_function => accel_function
    pc%outside_check => outside_check

    call pc%initialize(UC_elec_mass, 100*int(num_electrons_inception))
    call pc%use_cross_secs(max_energy, pdsim_table_size, cross_secs)

    where (pc%colls(:)%type == CS_ionize_t)
       pc%coll_is_event(:) = .true.
    end where

  end subroutine particles_initialize

  !> Check whether a particle is outside the gas
  integer function outside_check(my_part)
    type(PC_part_t), intent(inout) :: my_part
    integer                        :: i_cell, material
    real(dp)                       :: x(3)

    outside_check = 0

    if (pdsim_icdata_material > 0) then
       x = pdsim_convert_r(my_part%x)

       i_cell = 0
       material = -1
       call iu_get_icell_scalar_at(pdsim_ug, x, pdsim_icdata_material, &
            i_cell, material)

       if (material /= pdsim_gas_material_value) then
          outside_check = 1
       end if
    end if

  end function outside_check

  !> Get acceleration of particle by interpolating the electric field
  function accel_function(my_part) result(a)
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: a(3), x(3)

    x = pdsim_convert_r(my_part%x)
    a(pdsim_ndim+1:) = 0.0_dp

    call iu_interpolate_at(pdsim_ug, x, pdsim_ndim, &
         pdsim_pdata_field(1:pdsim_ndim), a(1:pdsim_ndim), my_part%id)
    a = a * UC_elec_q_over_m

    if (my_part%id < 1) error stop "accel_function: interpolation error"
  end function accel_function

  !> After updating the particles, events such as ionization are stored. Here
  !> we compute photoionization based on these events.
  subroutine handle_events(pc, avalanche_size)
    type(PC_t), intent(inout) :: pc
    integer, intent(inout)    :: avalanche_size

    integer                     :: n, n_photons
    real(dp), allocatable, save :: coords(:, :)
    type(PC_part_t)             :: my_part
    real(dp), parameter         :: array_incr_fac = 2.0_dp

    n = max(pc%n_events, pc%n_part)

    if (.not. allocated(coords)) then
       n = nint(n * array_incr_fac)
       allocate(coords(3, n))
    else if (size(coords, 2) < n) then
       n = nint(n * array_incr_fac)
       deallocate(coords)
       allocate(coords(3, n))
    end if

    avalanche_size = avalanche_size + &
         count(pc%event_list(1:pc%n_events)%ctype == CS_ionize_t)

    if (photoi_enabled) then
       call photoi_from_events(pc%n_events, pc%event_list, pc%rng, &
            coords, n_photons)

       do n = 1, n_photons
          my_part%x = pdsim_convert_r(coords(:, n))

          if (outside_check(my_part) == 0) then
             my_part%v(:)   = 0.0_dp
             my_part%w      = 1.0_dp
             my_part%t_left = 0.0_dp
             my_part%a(:)   = 0.0_dp ! Will be set later
             call pc%add_part(my_part)
          end if
       end do
    end if

    pc%n_events = 0

  end subroutine handle_events

  !> Write the particle coordinates to a file that can be loaded in Visit
  subroutine particles_write(i_run)
    integer, intent(in)       :: i_run
    character(len=200)        :: fname
    integer                   :: my_unit, n

    write(fname, "(A,I0.5,A,I0.5,A)") trim(pdsim_output_name) // '_run', &
         i_run, ".3D"

    open(newunit=my_unit, file=trim(fname), action='write')
    write(my_unit, *) "X Y Z eV"
    do n = 1, pc%n_part
       write(my_unit, *) pc%particles(n)%x, &
            PC_v_to_en(pc%particles(n)%v, pc%mass) / UC_elec_volt
    end do
    close(my_unit)

    print *, "Wrote ", trim(fname)
  end subroutine particles_write

  subroutine particles_simulate(cfg)
    type(CFG_t), intent(inout) :: cfg

    integer, allocatable :: avalanche_size(:)
    real(dp)             :: r_start(3), vars(2), initial_energy
    integer              :: n_runs, i_cell

    call particles_initialize(cfg)
    call CFG_get(cfg, "particles%initial_position", r_start)
    call CFG_get(cfg, "particles%n_runs", n_runs)
    call CFG_get(cfg, "particles%initial_energy", initial_energy)

    allocate(avalanche_size(n_runs))

    i_cell = 0
    call iu_interpolate_at(pdsim_ug, r_start, 2, [i_p_m1, i_kstar], vars, i_cell)

    call particles_simulate_avalanches(n_runs, r_start, initial_energy, &
         avalanche_size)
    call write_avalanche_size(n_runs, r_start, avalanche_size, vars(1), vars(2))

  end subroutine particles_simulate

  !> Simulate electron avalanches (possibly with photoionization) starting
  !> from initial electrons at given locations
  subroutine particles_simulate_avalanches(n_runs, r_start, initial_energy, &
       avalanche_size)
    integer, intent(in)  :: n_runs !< Number of runs
    real(dp), intent(in) :: r_start(3) !< Start location
    real(dp), intent(in) :: initial_energy !< Energy in eV
    integer, intent(out) :: avalanche_size(n_runs) !< Total number of ionizations
    type(PC_part_t)      :: initial_electron
    integer              :: n_output, i_run
    integer              :: max_iterations = huge(1)
    integer              :: iteration
    real(dp)             :: dt, time, speed
    logical              :: output_particles_now

    ! Set properties for initial electron
    initial_electron%v(:)   = 0.0_dp
    initial_electron%w      = 1.0_dp
    initial_electron%t_left = 0.0_dp
    initial_electron%a(:)   = 0.0_dp
    initial_electron%x      = r_start

    speed = sqrt(2/3.0_dp * initial_energy * UC_elec_volt / UC_elec_mass)

    do i_run = 1, n_runs
       call pc%remove_particles()

       ! Sample random direction as a Maxwellian
       initial_electron%v(:) = speed * [pc%rng%normal(), &
            pc%rng%normal(), pc%rng%normal()]

       call pc%add_part(initial_electron)
       avalanche_size(i_run) = 1

       time     = 0.0_dp
       dt       = max_dt
       n_output = 0

       do iteration = 1, max_iterations
          ! Set acceleration
          call pc%set_accel()

          output_particles_now = (particles_output .and. &
               time >= n_output * particles_output_dt)

          if (output_particles_now) then
             n_output = n_output + 1
             call particles_write(n_output)
          end if

          call pc%advance_openmp(dt)
          call pc%after_mover(dt)
          call pc%clean_up()
          time = time + dt

          call handle_events(pc, avalanche_size(i_run))

          if (pc%n_part == 0) exit
          if (pc%n_part >= num_electrons_inception) exit
       end do

       if (modulo(i_run, ceiling(1e-2_dp * n_runs)) == 0 &
            .and. pdsim_verbosity > 0) then
          write(*, "(F6.1,A)") (i_run*1e2_dp)/n_runs, "%"
       end if
    end do

  end subroutine particles_simulate_avalanches

  subroutine write_avalanche_size(n_runs, r_start, avalanche_size, &
       p_m1, K_star)
    integer, intent(in)  :: n_runs
    real(dp), intent(in) :: r_start(3)
    integer, intent(in)  :: avalanche_size(n_runs)
    real(dp), intent(in) :: p_m1
    real(dp), intent(in) :: K_star

    character(len=200) :: fname
    integer            :: my_unit, i_run
    real(dp)           :: r(3)

    write(fname, "(A,I0.5,A,I0.5,A)") trim(pdsim_output_name) // &
         "_avalanche_size.txt"

    open(newunit=my_unit, file=trim(fname), action="write")
    write(my_unit, *) p_m1, K_star, r_start

    r(:) = 0.0_dp
    do i_run = 1, n_runs
       write(my_unit, *) avalanche_size(i_run)
    end do

    close(my_unit)

    print *, "Wrote ", trim(fname)
  end subroutine write_avalanche_size

end module m_particles
