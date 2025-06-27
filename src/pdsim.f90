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
  use m_pq

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

  call simulate_avalanches(cfg)

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
    integer                    :: n
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
    call iu_add_point_data(pdsim_ug, "avalanche_p0", i_p0)
    call iu_add_point_data(pdsim_ug, "avalanche_w", i_w)
    call iu_add_point_data(pdsim_ug, "avalanche_time", i_travel_time)
    call iu_add_point_data(pdsim_ug, "avalanche_x1", i_x1)
    call iu_add_point_data(pdsim_ug, "avalanche_x2", i_x2)
    call iu_add_point_data(pdsim_ug, "avalanche_x3", i_x3)

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
       pdsim_ug%point_data(n, i_travel_time) = y(pdsim_ndim+2, n_steps)

       ! Store final position
       pdsim_ug%point_data(n, i_x1) = y(1, n_steps)
       pdsim_ug%point_data(n, i_x2) = y(2, n_steps)
       if (pdsim_ndim == 3) then
          pdsim_ug%point_data(n, i_x3) = y(3, n_steps)
       else
          pdsim_ug%point_data(n, i_x3) = 0.0_dp
       end if

       ! Store alpha and eta
       field = pdsim_ug%point_data(n, pdsim_pdata_field(1:pdsim_ndim))
       td = LT_get_mcol(pd%lkptbl, norm2(field))
       pdsim_ug%point_data(n, i_alpha) = td(pdsim_col_alpha)
       pdsim_ug%point_data(n, i_eta) = td(pdsim_col_eta)

       ! Compute w according to Kendall's 1948 paper
       call compute_kendall_w(pdsim_ndim, nvar, n_steps, y(:, 1:n_steps), &
            y_field(:, 1:n_steps), w)
       pdsim_ug%point_data(n, i_w) = w
       pdsim_ug%point_data(n, i_p0) = 1 - exp(y(pdsim_ndim+1, n_steps))/w
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
  subroutine compute_kendall_w(ndim, nvar, n_steps, y, y_field, w)
    integer, intent(in)   :: ndim
    integer, intent(in)   :: nvar
    integer, intent(in)   :: n_steps
    real(dp), intent(in)  :: y(ndim+nvar, n_steps)
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

  !> Simulate avalanches as discrete events
  subroutine simulate_avalanches(cfg)
    use m_random
    type(cfg_t), intent(inout) :: cfg

    integer                        :: n, ix, k, i_run, i_avalanche
    integer                        :: max_avalanches, n_runs
    integer                        :: max_photons
    integer                        :: inception_threshold
    integer                        :: n_photons, i_cell
    integer, parameter             :: n_vars = 6
    integer                        :: i_vars(n_vars)
    real(dp)                       :: r(3), w, k_integral, travel_time
    real(dp)                       :: time, r_arrival(3)
    real(dp)                       :: vars(n_vars)
    real(dp)                       :: p_avg, volume
    type(avalanche_t), allocatable :: avalanches(:)
    real(dp), allocatable          :: absorption_locations(:, :)
    real(dp), allocatable          :: inception_time(:)
    logical, allocatable           :: inception(:)
    type(pqr_t)                    :: pq
    type(rng_t)                    :: rng

    i_vars = [i_w, i_k_integral, i_travel_time, i_x1, i_x2, i_x3]

    call iu_add_point_data(pdsim_ug, "inception_time", i_inception_time)
    call iu_add_point_data(pdsim_ug, "inception_prob", i_inception_prob)

    call CFG_get(cfg, "avalanche%max_total_number", max_avalanches)
    call CFG_get(cfg, "avalanche%max_photons", max_photons)
    call CFG_get(cfg, "avalanche%n_runs", n_runs)
    call CFG_get(cfg, "avalanche%inception_threshold", inception_threshold)

    allocate(avalanches(max_avalanches))
    allocate(absorption_locations(3, max_photons))
    allocate(inception(n_runs))
    allocate(inception_time(n_runs))

    ! Create priority queue that will store upcoming avalanches, sorted by
    ! their arrival time
    call pqr_create(pq)

    call rng%set_random_seed()

    do n = 1, pdsim_ug%n_points
       if (modulo(n, pdsim_ug%n_points/100) == 0) then
          write(*, "(F6.1,A)") (n*1e2_dp)/pdsim_ug%n_points, "%"
       end if

       do i_run = 1, n_runs
          i_avalanche = 0
          time = 0.0_dp
          pq%n_stored = 0

          ! Parameters for the initial avalanche
          r = pdsim_ug%points(:, n)
          w = pdsim_ug%point_data(n, i_w)
          k_integral = pdsim_ug%point_data(n, i_k_integral)
          travel_time = pdsim_ug%point_data(n, i_travel_time)
          r_arrival = pdsim_ug%point_data(n, [i_x1, i_x2, i_x3])

          call add_new_avalanche(rng, time, r, w, k_integral, travel_time, &
               r_arrival, i_avalanche, avalanches, pq)

          time_loop: do while (pq%n_stored > 0)
             ! Get the next avalanche
             call pqr_pop(pq, ix, time)

             associate (av => avalanches(ix))
               ! Produce photons
               call photoi_sample_photons(rng, av%r_source, &
                    real(av%avalanche_size, dp), n_photons, absorption_locations)

               do k = 1, n_photons
                  i_cell = 0
                  r = pdsim_convert_r(absorption_locations(:, k))
                  call iu_get_cell(pdsim_ug, r, i_cell)

                  if (i_cell > 0) then
                     ! Photoelectron can contribute to new avalanche
                     call iu_interpolate_at(pdsim_ug, r, n_vars, &
                          i_vars, vars, i_cell)

                     w = vars(1)
                     k_integral = vars(2)
                     travel_time = vars(3)
                     r_arrival = vars(4:6)

                     call add_new_avalanche(rng, time, r, w, k_integral, &
                          travel_time, r_arrival, i_avalanche, avalanches, pq)

                     ! Exit when the inception threshold has been reached
                     if (pq%n_stored == inception_threshold) exit time_loop
                  end if
               end do
             end associate
          end do time_loop

          inception(i_run) = (pq%n_stored >= inception_threshold)
          inception_time(i_run) = time
       end do

       pdsim_ug%point_data(n, i_inception_time) = &
            sum(inception_time, mask=inception)/max(1, count(inception))
       pdsim_ug%point_data(n, i_inception_prob) = count(inception) / real(n_runs, dp)
    end do

    call compute_pointdata_average(pdsim_ug, i_inception_prob, &
         pdsim_axisymmetric, p_avg, volume)

    write(*, "(A,E11.3)") " Average inception probability: ", p_avg
    write(*, "(A,E12.4)") " Total volume of gas: ", volume

  end subroutine simulate_avalanches

  !> Add a new avalanche to the list of avalanches, if it has a nonzero size
  subroutine add_new_avalanche(rng, time, r, w, k_integral, &
       dt, r_arrival, ix, avalanches, pq)
    use m_random
    use iso_fortran_env, only: int64
    type(rng_t), intent(inout)       :: rng
    real(dp), intent(in)             :: time
    real(dp), intent(in)             :: r(3)
    real(dp), intent(in)             :: w
    real(dp), intent(in)             :: k_integral
    real(dp), intent(in)             :: dt
    real(dp), intent(in)             :: r_arrival(3)
    integer, intent(inout)           :: ix
    type(avalanche_t), intent(inout) :: avalanches(:)
    type(pqr_t), intent(inout)       :: pq

    real(dp) :: p0, pgeom, tmp

    ! Probability of the avalanche having size zero
    p0 = 1 - exp(k_integral)/w

    if (rng%unif_01() > p0) then
       ! Add avalanche
       ix = ix + 1

       if (ix > size(avalanches)) &
            error stop "Not enough storage for avalanches"

       ! Probability of geometric distribution
       pgeom = 1/w

       associate (av => avalanches(ix))
         av%t_source = time
         av%r_source = r
         av%r_arrival = r_arrival
         av%t_arrival = time + dt

         ! Sample avalanche size. Take care of cases when pgeom >= 1.0 due to
         ! numerical errors
         if (pgeom < 1) then
            tmp = log(1 - rng%unif_01()) / log(1 - pgeom)
         else
            tmp = 0.0_dp
         end if

         av%avalanche_size = ceiling(tmp, int64)

         call pqr_push(pq, ix, av%t_arrival)
       end associate
    end if

  end subroutine add_new_avalanche

  !> Compute volume average of a variable defined at points (vertices)
  subroutine compute_pointdata_average(ug, iv, axisymmetric, avg, vol)
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

  end subroutine compute_pointdata_average

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

    call write_ionization_count(n_pos, r_start, n_ionizations)

  end subroutine particle_simulation

  subroutine write_ionization_count(n_pos, r_start, n_ionizations)
    integer, intent(in)  :: n_pos
    real(dp), intent(in) :: r_start(pdsim_ndim, n_pos)
    integer, intent(in)  :: n_ionizations(n_pos)

    character(len=200) :: fname
    integer            :: my_unit, n
    real(dp)           :: r(3)

    write(fname, "(A,I0.5,A,I0.5,A)") trim(pd%output_name) // &
         "_ionization_count.3D"

    open(newunit=my_unit, file=trim(fname), action="write")
    write(my_unit, *) "X Y Z n_ionizations"

    r(:) = 0.0_dp
    do n = 1, n_pos
       r(1:pdsim_ndim) = r_start(:, n)
       write(my_unit, *) r, n_ionizations(n)
    end do

    close(my_unit)

    print *, "Wrote ", trim(fname)
  end subroutine write_ionization_count

end program simulate_inception
