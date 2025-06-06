program simulate_inception
  use m_config
  use m_pdsim
  use m_cross_sec
  use m_particle_core
  use m_units_constants
  use m_gas
  use m_photoi
  use m_lookup_table
  use m_interp_unstructured

  implicit none

  type(CFG_t)   :: cfg
  type(pdsim_t) :: pd

  integer, allocatable  :: i_start_cells(:)
  real(dp), allocatable :: r_start(:, :), tmp_array(:)
  integer, allocatable  :: n_ionizations(:)
  integer               :: i, n_pos, var_size

  call pdsim_create_config(cfg)

  call CFG_update_from_arguments(cfg)
  call CFG_check(cfg)

  call pdsim_initialize(cfg, pd)

  if (pd%compute_ionization_integral) then
     call compute_ionization_integral(cfg)
  end if

  select case (pd%start_method)
  case ("random_cell")
     if (pd%n_random_positions == -1) &
          error stop "simulation%n_random_positions not specified"
     call pdsim_get_gas_cells(pdsim_ug, pd%n_random_positions, &
          pd%initial_rmin, pd%initial_rmax, i_start_cells)
     allocate(r_start(pdsim_ndim, pd%n_random_positions))
     do i = 1, pd%n_random_positions
        r_start(:, i) = iu_get_cell_center(pdsim_ug, i_start_cells(i))
     end do
     n_pos = pd%n_random_positions
  case ("location")
     call CFG_get_size(cfg, "simulation%initial_positions", var_size)
     if (modulo(var_size, pdsim_ndim) /= 0) &
          error stop "simulation%initial_positions has invalid size"
     if (var_size == 0) error stop "simulation%initial_positions not specified"
     n_pos = var_size/pdsim_ndim

     allocate(tmp_array(var_size))
     call CFG_get(cfg, "simulation%initial_positions", tmp_array)
     r_start = reshape(tmp_array, [pdsim_ndim, n_pos])
  case default
     error stop "simulation%start_method can be random_cell, location"
  end select


  if (pd%simulate_particles) then
     allocate(n_ionizations(n_pos))
     n_ionizations(:) = 0
     call particle_simulation(n_pos, r_start, n_ionizations)
  end if

  call iu_write_vtk(pdsim_ug, trim(pd%output_name) // ".vtu")
  print *, "Wrote ", trim(pd%output_name) // ".vtu"

contains

  subroutine compute_ionization_integral(cfg)
    type(cfg_t), intent(inout) :: cfg
    integer                    :: n, i_k_integral, i_alpha, i_eta
    integer                    :: max_steps, n_steps
    real(dp)                   :: min_dx, max_dx, boundary_distance
    real(dp)                   :: rtol, atol, r(3), y(pdsim_ndim+1)
    real(dp)                   :: field(pdsim_ndim), alpha_eta(2)
    logical                    :: reverse

    call CFG_get(cfg, "integral%max_steps", max_steps)
    call CFG_get(cfg, "integral%rtol", rtol)
    call CFG_get(cfg, "integral%atol", atol)
    call CFG_get(cfg, "integral%min_dx", min_dx)
    call CFG_get(cfg, "integral%max_dx", max_dx)
    call CFG_get(cfg, "integral%boundary_distance", boundary_distance)
    reverse = .true.

    call iu_add_point_data(pdsim_ug, "K_integral", i_k_integral)
    call iu_add_point_data(pdsim_ug, "alpha", i_alpha)
    call iu_add_point_data(pdsim_ug, "eta", i_eta)

    !$omp parallel do private(r, y, n_steps, field, alpha_eta)
    do n = 1, pdsim_ug%n_points
       r = pdsim_ug%points(:, n)

       ! Ensure points are some distance away from boundaries
       r(1:pdsim_ndim) = max(r(1:pdsim_ndim), &
            pdsim_ug%rmin(1:pdsim_ndim) + boundary_distance)
       r(1:pdsim_ndim) = min(r(1:pdsim_ndim), &
            pdsim_ug%rmax(1:pdsim_ndim) - boundary_distance)

       call iu_integrate_along_field(pdsim_ug, alpha_eff_func, pdsim_ndim, &
            r(1:pdsim_ndim), pdsim_pdata_field(1:pdsim_ndim), &
            min_dx, max_dx, max_steps, rtol, atol, reverse, y, n_steps, &
            pdsim_cdata_material, pdsim_gas_material_value)

       pdsim_ug%point_data(n, i_k_integral) = y(3)

       ! Store alpha and eta
       field = pdsim_ug%point_data(n, pdsim_pdata_field(1:pdsim_ndim))
       alpha_eta = LT_get_mcol(pd%lkptbl, norm2(field))
       pdsim_ug%point_data(n, i_alpha) = alpha_eta(1)
       pdsim_ug%point_data(n, i_eta) = alpha_eta(2)
    end do
    !$omp end parallel do

  end subroutine compute_ionization_integral

  function alpha_eff_func(ndim, r, field) result(alpha_eff)
    integer, intent(in)  :: ndim
    real(dp), intent(in) :: r(ndim)
    real(dp), intent(in) :: field(ndim)
    real(dp)             :: alpha_eff, alpha_eta(2)

    alpha_eta = LT_get_mcol(pd%lkptbl, norm2(field))
    alpha_eff = alpha_eta(1) - alpha_eta(2)
  end function alpha_eff_func

  subroutine particle_simulation(n_pos, r_start, n_ionizations)
    integer, intent(in)  :: n_pos
    real(dp), intent(in) :: r_start(pdsim_ndim, n_pos)
    integer, intent(out) :: n_ionizations(n_pos)
    type(PC_part_t)      :: initial_electron
    integer              :: n_output, i_start
    integer              :: max_iterations = huge(1)
    integer              :: iteration
    real(dp)             :: dt, time
    logical              :: output_particles_now

    ! Set properties for initial electron
    initial_electron%v(:)   = 0.0_dp
    initial_electron%w      = 1.0_dp
    initial_electron%t_left = 0.0_dp
    initial_electron%a(:)   = 0.0_dp
    initial_electron%x(:)   = 0.0_dp

    do i_start = 1, n_pos
       call pd%pc%remove_particles()

       initial_electron%x(1:pdsim_ndim) = r_start(:, i_start)
       call pd%pc%add_part(initial_electron)

       time     = 0.0_dp
       dt       = pd%max_dt
       n_output = 0

       do iteration = 1, max_iterations
          ! Set acceleration
          call pd%pc%set_accel()

          output_particles_now = (pd%output_particles .and. &
               time >= n_output * pd%output_particles_dt)

          if (output_particles_now) then
             n_output = n_output + 1
             call pdsim_write_particles(pd, n_output, i_start)
          end if

          call pd%pc%advance(dt)
          call pd%pc%after_mover(dt)
          call pd%pc%clean_up()
          time = time + dt

          call pdsim_handle_events(pd%pc, n_ionizations(i_start))

          if (pd%pc%n_part == 0) exit
          if (pd%pc%n_part >= pd%num_electrons_inception) exit
       end do

       print *, i_start, iteration, n_ionizations(i_start)
    end do

    call write_ionization_count(i_start_cells, n_ionizations)

  end subroutine particle_simulation

  subroutine write_ionization_count(i_start_cells, n_ionizations)
    integer, intent(in) :: i_start_cells(:)
    integer, intent(in) :: n_ionizations(:)

    character(len=200)        :: fname
    integer                   :: my_unit, n

    write(fname, "(A,I0.5,A,I0.5,A)") trim(pd%output_name) // &
         "_ionization_count.3D"

    open(newunit=my_unit, file=trim(fname), action="write")
    write(my_unit, *) "X Y Z n_ionizations"
    do n = 1, size(i_start_cells)
       write(my_unit, *) iu_get_cell_center(pdsim_ug, i_start_cells(n)), &
            n_ionizations(n)
    end do
    close(my_unit)

    print *, "Wrote ", trim(fname)
  end subroutine write_ionization_count

end program simulate_inception
