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
    integer                    :: i_pgeom, i_x1, i_x2, i_x3, i_time
    integer                    :: max_steps, n_steps, nvar
    real(dp)                   :: min_dx, max_dx, boundary_distance
    real(dp)                   :: rtol, atol, r(3), w
    real(dp)                   :: field(pdsim_ndim), td(pdsim_ncols)
    real(dp)                   :: domain_center(3), direction(3)
    real(dp), allocatable      :: y(:, :), y_field(:, :)
    logical                    :: reverse

    call CFG_get(cfg, "integral%max_steps", max_steps)
    call CFG_get(cfg, "integral%rtol", rtol)
    call CFG_get(cfg, "integral%atol", atol)
    call CFG_get(cfg, "integral%min_dx", min_dx)
    call CFG_get(cfg, "integral%max_dx", max_dx)
    call CFG_get(cfg, "integral%boundary_distance", boundary_distance)
    reverse = .true.
    nvar = 2

    call iu_add_point_data(pdsim_ug, "K_integral", i_k_integral)
    call iu_add_point_data(pdsim_ug, "alpha", i_alpha)
    call iu_add_point_data(pdsim_ug, "eta", i_eta)
    call iu_add_point_data(pdsim_ug, "avalanche_pgeom", i_pgeom)
    call iu_add_point_data(pdsim_ug, "avalanche_time", i_time)
    call iu_add_point_data(pdsim_ug, "avalanche_x1", i_x1)
    call iu_add_point_data(pdsim_ug, "avalanche_x2", i_x2)

    if (pdsim_ndim == 3) then
       call iu_add_point_data(pdsim_ug, "avalanche_x3", i_x3)
    end if

    domain_center = 0.5_dp * (pdsim_ug%rmin + pdsim_ug%rmax)

    allocate(y(pdsim_ndim+nvar, max_steps))
    allocate(y_field(pdsim_ndim, max_steps))

    !$omp parallel do private(r, y, y_field, n_steps, field, td, &
    !$omp w, direction)
    do n = 1, pdsim_ug%n_points
       r = pdsim_ug%points(:, n)

       ! Ensure points are some distance away from boundaries
       if (pdsim_ug%point_is_at_boundary(n)) then
          direction = (domain_center - r) / norm2(domain_center - r)
          r = r + boundary_distance * direction
       end if

       call iu_integrate_along_field(pdsim_ug, pdsim_ndim, alpha_eff_sub, &
            r(1:pdsim_ndim), pdsim_pdata_field(1:pdsim_ndim), &
            min_dx, max_dx, max_steps, rtol, atol, reverse, nvar, y, y_field, &
            n_steps, pdsim_axisymmetric, &
            pdsim_cdata_material, pdsim_gas_material_value)

       if (n_steps > max_steps) then
          print *, "Error for start position", r(1:pdsim_ndim)
          print *, "n_steps = ", n_steps, ", max_steps = ", max_steps
          error stop "Increase integral%max_steps"
       end if

       pdsim_ug%point_data(n, i_k_integral) = y(pdsim_ndim+1, n_steps)
       pdsim_ug%point_data(n, i_time) = y(pdsim_ndim+2, n_steps)

       ! Store final position
       pdsim_ug%point_data(n, i_x1) = y(1, n_steps)
       pdsim_ug%point_data(n, i_x2) = y(2, n_steps)
       if (pdsim_ndim == 3) then
          pdsim_ug%point_data(n, i_x3) = y(3, n_steps)
       end if

       ! Store alpha and eta
       field = pdsim_ug%point_data(n, pdsim_pdata_field(1:pdsim_ndim))
       td = LT_get_mcol(pd%lkptbl, norm2(field))
       pdsim_ug%point_data(n, i_alpha) = td(pdsim_col_alpha)
       pdsim_ug%point_data(n, i_eta) = td(pdsim_col_eta)

       ! Compute w according to Kendall's 1948 paper
       call compute_kendall_w(pdsim_ndim, n_steps, y(:, 1:n_steps), &
            y_field(:, 1:n_steps), w)
       pdsim_ug%point_data(n, i_pgeom) = 1/w
    end do
    !$omp end parallel do

  end subroutine compute_ionization_integral

  !> Computes alpha effective and 1/velocity
  subroutine alpha_eff_sub(ndim, r, field, nvar, integrand)
    integer, intent(in)   :: ndim
    real(dp), intent(in)  :: r(ndim)
    real(dp), intent(in)  :: field(ndim)
    integer, intent(in)   :: nvar
    real(dp), intent(out) :: integrand(nvar)
    real(dp)              :: field_norm, td(pdsim_ncols)

    field_norm = norm2(field)
    td = LT_get_mcol(pd%lkptbl, field_norm)
    integrand(1) = td(pdsim_col_alpha) - td(pdsim_col_eta)
    integrand(2) = 1/(field_norm * td(pdsim_col_mu))
  end subroutine alpha_eff_sub

  !> Compute the w value at the final point, as defined in Kendall's 1948
  !> paper. This can be written as:
  !> w = 1 + exp(-rho(x)) * integral(exp(rho(x')) * alpha(x') dx')
  !> for x' between 0 and x. Here rho(x) = -K(x).
  subroutine compute_kendall_w(ndim, n_steps, y, y_field, w)
    integer, intent(in)   :: ndim
    integer, intent(in)   :: n_steps
    real(dp), intent(in)  :: y(ndim+1, n_steps)
    real(dp), intent(in)  :: y_field(ndim, n_steps)
    real(dp), intent(out) :: w

    integer               :: n
    real(dp)              :: alpha, dx, rho
    real(dp)              :: integrand(n_steps)

    ! Determine integrand at every point
    do n = 1, n_steps
       alpha = LT_get_col(pd%lkptbl, pdsim_col_alpha, norm2(y_field(:, n)))
       rho = -y(ndim+1, n)
       integrand(n) = exp(rho) * alpha
    end do

    ! Use composite trapezoidal rule to evaluate integral
    w = 0.0_dp
    do n = 1, n_steps-1
       dx = norm2(y(1:ndim, n+1) - y(1:ndim, n))
       w = w + dx * 0.5_dp * (integrand(n) + integrand(n+1))
    end do

    w = 1 + exp(y(ndim+1, n_steps)) * w
  end subroutine compute_kendall_w

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
