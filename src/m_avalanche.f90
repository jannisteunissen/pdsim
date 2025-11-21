!> Module for simulation of avalanches
module m_avalanche
  use iso_fortran_env, only: int64
  use m_config
  use m_interp_unstructured
  use m_pdsim
  use m_photoi
  use m_pq
  use m_random

  implicit none
  private

  type avalanche_t
     !> Start time
     real(dp)       :: t_source
     !> Start location
     real(dp)       :: r_source(3)
     !> Arrival time
     real(dp)       :: t_arrival
     !> Arrival location of avalanche
     real(dp)       :: r_arrival(3)
     !> Arrival location of positive ions
     real(dp)       :: r_ion_arrival(3)
     !> Secondary emission coefficient for ions
     real(dp)       :: ion_gamma
     !> Ion travel time for secondary emission
     real(dp)       :: ion_travel_time
     !> Number of ionizations produced by the initial electron
     integer(int64) :: num_ionizations
  end type avalanche_t

  integer, parameter :: n_vars = 11
  integer            :: i_vars(n_vars)

  !> Produce photoemission this distance away from the boundary
  real(dp) :: photoemission_boundary_distance

  ! Public methods
  public :: avalanche_create_config
  public :: avalanche_simulate

contains

  subroutine avalanche_create_config(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "avalanche%n_runs", 10, &
         "Number of runs per initial avalanche location")
    call CFG_add(cfg, "avalanche%use_antithetic", .true., &
         "Whether to use antithetic variates to reduce variance")
    call CFG_add(cfg, "avalanche%use_early_exit", .false., &
         "Exit after finding point with non-zero inception probability")
    call CFG_add(cfg, "avalanche%max_photons", 100*1000, &
         "Maximum number of photons produced by a single avalanche")
    call CFG_add(cfg, "avalanche%inception_count", 1000, &
         "Inception occurs when there are this many future avalanches")
    call CFG_add(cfg, "avalanche%inception_size", 1e8_dp, &
         "Inception occurs when an avalanche of at least this size occurs")
    call CFG_add(cfg, "avalanche%trace_photons", .false., &
         "Whether to trace photon paths")
    call CFG_add(cfg, "avalanche%limit_photoelectrons", .true., &
         "Limit number of photoelectrons produced by a single avalanche")
    call CFG_add(cfg, "avalanche%photoemission_boundary_distance", 1e-9_dp, &
         "Produce photoemission this distance away from the boundary (m)")
    call CFG_add(cfg, "avalanche%r_min", [-1e100_dp, -1e100_dp, -1e100_dp], &
         "Only start avalanches in a box between r_min and r_max")
    call CFG_add(cfg, "avalanche%r_max", [1e100_dp, 1e100_dp, 1e100_dp], &
         "Only start avalanches in a box between r_min and r_max")

  end subroutine avalanche_create_config

  !> Simulate avalanches as discrete events
  subroutine avalanche_simulate(cfg)
    use omp_lib
    type(cfg_t), intent(inout) :: cfg

    integer                        :: n
    integer                        :: n_runs, thread_id
    integer                        :: max_photons
    integer                        :: inception_count
    integer(int64)                 :: n_steps_total, n_steps_point
    real(dp)                       :: inception_size
    real(dp)                       :: p_avg, pvar_avg, volume
    real(dp)                       :: r_min(3), r_max(3)
    logical                        :: trace_photons, use_antithetic
    logical                        :: limit_photoelectrons
    logical                        :: use_early_exit
    type(avalanche_t), allocatable :: avalanches(:)
    logical, allocatable           :: skip_point(:)
    real(dp)                       :: inception_time
    real(dp)                       :: inception_probability
    real(dp)                       :: inception_pvar
    type(pqr_t)                    :: pq
    type(rng_t)                    :: rng
    type(prng_t)                   :: prng
    integer(int64)                 :: t_start, t_end, count_rate

    call iu_reserve_point_data_storage(pdsim_ug, 3)
    call iu_add_point_data(pdsim_ug, "inception_time", i_inception_time)
    call iu_add_point_data(pdsim_ug, "inception_prob", i_inception_prob)
    call iu_add_point_data(pdsim_ug, "inception_pvar", i_inception_pvar)

    call CFG_get(cfg, "avalanche%max_photons", max_photons)
    call CFG_get(cfg, "avalanche%n_runs", n_runs)
    call CFG_get(cfg, "avalanche%use_antithetic", use_antithetic)
    call CFG_get(cfg, "avalanche%use_early_exit", use_early_exit)
    call CFG_get(cfg, "avalanche%inception_count", inception_count)
    call CFG_get(cfg, "avalanche%inception_size", inception_size)
    call CFG_get(cfg, "avalanche%trace_photons", trace_photons)
    call CFG_get(cfg, "avalanche%limit_photoelectrons", limit_photoelectrons)
    call CFG_get(cfg, "avalanche%photoemission_boundary_distance", &
         photoemission_boundary_distance)
    call CFG_get(cfg, "avalanche%r_min", r_min)
    call CFG_get(cfg, "avalanche%r_max", r_max)

    if (use_early_exit) then
       if (.not. omp_get_cancellation()) then
          error stop "Set environment variable OMP_CANCELLATION to true"
       end if
    end if

    allocate(avalanches(inception_count))
    allocate(skip_point(pdsim_ug%n_points))

    do n = 1, pdsim_ug%n_points
       skip_point(n) = any(pdsim_ug%points(:, n) < r_min .or. &
            pdsim_ug%points(:, n) > r_max)
    end do

    ! Set module-level variable
    i_vars = [i_p_m1, i_kstar, i_avalanche_time, i_x1, i_x2, i_x3, &
         i_ion_x1, i_ion_x2, i_ion_x3, i_ion_gamma, i_ion_time]

    ! Create priority queue that will store upcoming avalanches, sorted by
    ! their arrival time
    call pqr_create(pq, inception_count)

    call rng%set_random_seed()
    call prng%init_parallel(omp_get_max_threads(), rng)

    !$omp parallel private(n, thread_id, pq, avalanches, n_steps_point, &
    !$omp inception_probability, inception_time, inception_pvar)

    ! Get a random number generator for each thread
    thread_id = omp_get_thread_num() + 1
    rng = prng%rngs(thread_id)
    n_steps_total = 0

    call system_clock(t_start, count_rate)
    !$omp do schedule(dynamic) reduction(+:n_steps_total) 
    do n = 1, pdsim_ug%n_points
       if (modulo(n, ceiling(5e-2_dp * pdsim_ug%n_points)) == 0 .and. &
            pdsim_verbosity > 0) then
          write(*, "(F6.1,A)") (n*1e2_dp)/pdsim_ug%n_points, "%"
       end if

       if (skip_point(n)) then
          pdsim_ug%point_data(n, i_inception_time) = 0.0_dp
          pdsim_ug%point_data(n, i_inception_prob) = 0.0_dp
          pdsim_ug%point_data(n, i_inception_pvar) = 0.0_dp
          cycle
       end if

       call run_avalanches(n, n_runs, inception_count, inception_size, &
            max_photons, trace_photons, limit_photoelectrons, use_antithetic, &
            rng, pq, avalanches, inception_probability, inception_time, &
            inception_pvar, n_steps_point)

       n_steps_total = n_steps_total + n_steps_point
       pdsim_ug%point_data(n, i_inception_time) = inception_time
       pdsim_ug%point_data(n, i_inception_prob) = inception_probability
       pdsim_ug%point_data(n, i_inception_pvar) = inception_pvar

       if (use_early_exit .and. inception_probability > 0) then
          !$omp cancel do
       end if
    end do
    !$omp end do
    !$omp end parallel
    call system_clock(t_end, count_rate)

    call pdsim_pointdata_average(pdsim_ug, i_inception_prob, &
         pdsim_axisymmetric, 1, p_avg, volume)
    call pdsim_pointdata_average(pdsim_ug, i_inception_pvar, &
         pdsim_axisymmetric, 2, pvar_avg, volume)

    if (pdsim_verbosity > 0) then
       if (use_early_exit .and. p_avg > 0) then
          print *, "The simulation was stopped after detecting inception"
       else
          print *, "The simulation completed without detecting inception"
       end if

       write(*, "(A,E12.4)") " Average inception probability: ", p_avg
       write(*, "(A,E12.4)") " Standard deviation bound:      ", sqrt(pvar_avg)
       write(*, "(A,E12.4)") " Total volume of gas:           ", volume
       write(*, "(A,I0)")    " Total number of steps taken:     ", n_steps_total
       write(*, "(A,E12.4)") " Time for avalanches (s):       ", &
            (t_end - t_start)/real(count_rate, dp)
    else
       write(*, "(3E12.4)") p_avg, sqrt(pvar_avg), volume
    end if

  end subroutine avalanche_simulate

  !> Run n_runs avalanches starting at a point
  subroutine run_avalanches(ip, n_runs_max, inception_count, inception_size, &
       max_photons, trace_photons, limit_photoelectrons, use_antithetic, rng, &
       pq, avalanches, inception_probability, inception_time, inception_pvar, &
       n_steps_point)
    integer, intent(in)              :: ip !< Point index
    integer, intent(in)              :: n_runs_max
    integer, intent(in)              :: inception_count
    real(dp), intent(in)             :: inception_size
    integer, intent(in)              :: max_photons
    logical, intent(in)              :: trace_photons
    logical, intent(in)              :: limit_photoelectrons
    logical, intent(in)              :: use_antithetic
    type(rng_t), intent(inout)       :: rng
    type(pqr_t), intent(inout)       :: pq
    type(avalanche_t), intent(inout) :: avalanches(inception_count)
    real(dp), intent(out) :: inception_probability
    real(dp), intent(out) :: inception_time
    real(dp), intent(out) :: inception_pvar
    integer(int64), intent(out) :: n_steps_point

    logical             :: inception(n_runs_max)
    real(dp)            :: t_inception(n_runs_max)
    ! Factor for inception probability
    real(dp)            :: p_factor

    integer, parameter :: max_steps = 1000*1000
    integer           :: ix, n_runs, i_run, i_half, i_step
    real(dp)          :: time, r(3), vars(n_vars), p, p_bnd
    real(dp)          :: p_m1, rand_num(n_runs_max), num_expected
    type(avalanche_t) :: av

    inception(:) = .false.
    t_inception(:) = 0.0_dp
    n_steps_point = 0

    ! Parameters for the initial avalanche
    r = pdsim_ug%points(:, ip)
    vars = pdsim_ug%point_data(ip, i_vars)

    ! Determine expected number of avalanches producing additional ionization
    p_m1 = vars(1)
    num_expected = (1 - p_m1) * n_runs_max
    n_runs = min(ceiling(num_expected), n_runs_max)

    do i_run = 1, n_runs
       rand_num(i_run) = rng%unif_01()
    end do

    if (use_antithetic) then
       ! Use antithetic variables for second half of runs
       i_half = (n_runs+1)/2

       do i_run = i_half + 1, n_runs
          ! Avoid rare case of exactly 1.0_dp to avoid issues with a log(1-x)
          rand_num(i_run) = 1.0_dp - max(rand_num(i_run - i_half), &
               epsilon(1.0_dp))
       end do
    end if

    do i_run = 1, n_runs
       time = 0.0_dp
       call pqr_reset(pq)

       call add_new_avalanche(rng, time, r, vars, avalanches, pq, &
            u01_created=1.0_dp, u01_size=rand_num(i_run))

       do i_step = 1, max_steps

          if (pq%n_stored == 0) then
             exit
          else
             ! Get the next avalanche
             call pqr_pop_aix(pq, ix, time)
          end if

          ! Store copy of avalanche, as the index might be overwritten
          av = avalanches(ix)

          ! Check if avalanche size exceeds threshold
          if (avalanches(ix)%num_ionizations > inception_size) then
             inception(i_run) = .true.
             exit
          end if

          if (photoi_enabled) then
             call photoi_SEE(av, time, max_photons, trace_photons, &
                  limit_photoelectrons, inception_count, rng, avalanches, pq)

             if (pq%n_stored >= inception_count) then
                inception(i_run) = .true.
                exit
             end if
          end if

          if (pdsim_ion_see_enabled) then
             call ion_SEE(av, time, inception_count, &
                  rng, avalanches, pq)

             if (pq%n_stored >= inception_count) then
                inception(i_run) = .true.
                exit
             end if
          end if
       end do

       if (inception(i_run)) then
          ! Store inception time
          t_inception(i_run) = time
       end if

       n_steps_point = n_steps_point + i_step
    end do

    if (n_runs > 0) then
       ! Correction factor for inception probability
       p_factor = num_expected / n_runs_max

       ! p is an estimate of the probability per actual run
       p = count(inception) / real(n_runs, dp)

       ! Estimate upper bound for sample variance of p
       p_bnd = 0.5_dp
       inception_pvar = p_bnd * (1 - p_bnd) / n_runs * p_factor**2
    else
       p_factor = 0.0_dp
       p = 0.0_dp
       inception_pvar = 0.0_dp
    end if

    inception_time = sum(t_inception, mask=inception)/max(1, count(inception))
    inception_probability = p * p_factor

  end subroutine run_avalanches

  !> Sample photoionization secondary electron emission from an avalanche and
  !> store the resulting new avalanches
  subroutine photoi_SEE(av, time, max_photons, trace_photons, &
       limit_photoelectrons, inception_count, rng, avalanches, pq)
    type(avalanche_t), intent(in)    :: av
    real(dp), intent(in)             :: time
    integer, intent(in)              :: max_photons
    logical, intent(in)              :: trace_photons
    logical, intent(in)              :: limit_photoelectrons
    integer, intent(in)              :: inception_count
    type(rng_t), intent(inout)       :: rng
    type(avalanche_t), intent(inout) :: avalanches(:)
    type(pqr_t), intent(inout)       :: pq

    real(dp) :: absorption_locations(3, max_photons)
    integer  :: k, n, n_photons, i_cell, i_start, status, n_steps, material
    integer  :: n_stored_limit, n_photoemission
    real(dp) :: r(3), r_p(3), vars(n_vars), dvec(3), gamma

    ! Subtract one since initial ionization does not produce photons
    call photoi_sample_photons(rng, av%r_arrival, &
         av%num_ionizations - 1.0_dp, max_photons, &
         n_photons, absorption_locations)

    if (n_photons == 0) return

    if (limit_photoelectrons) then
       ! Limit for number of stored avalanches. Prevents a single avalanche
       ! from producing so many photons that the inception criterion is met.
       n_stored_limit = min(inception_count, pq%n_stored + inception_count/2)
    else
       n_stored_limit = inception_count
    end if

    if (trace_photons) then
       ! Get start location
       i_start = 0
       call iu_get_cell(pdsim_ug, av%r_arrival, i_start)

       ! Due to interpolation, av%r_arrival might not lie inside the gas. In
       ! such cases, do not create photons
       if (i_start == 0) then
          return
       else if (pdsim_ug%icell_data(i_start, pdsim_icdata_material) /= &
            pdsim_gas_material_value) then
          return
       end if

       do k = 1, n_photons
          r = pdsim_convert_r(absorption_locations(:, k))

          ! Trace path to destination
          call iu_get_cell_through_neighbors(pdsim_ug, huge(1), av%r_arrival, &
               r, i_start, i_cell, r_p, n_steps, status, pdsim_icdata_material)

          if (status == 0) then
             ! Photoelectron in gas can contribute to new avalanche
             call iu_interpolate_at(pdsim_ug, r, n_vars, &
                  i_vars, vars, i_cell)

             call add_new_avalanche(rng, time, r, vars, avalanches, pq)

          else if (pdsim_photoemission_enabled) then
             if (status == 1) then
                ! Photon hit another material
                material = pdsim_ug%icell_data(i_cell, pdsim_icdata_material)
                gamma = pdsim_photoemission_gamma(1+material)
             else
                ! Photon hit a domain boundary
                gamma = pdsim_photoemission_gamma(1)
             end if

             n_photoemission = rng%poisson(gamma)

             if (n_photoemission > 0) then
                ! Create electron slightly away from boundary
                dvec = av%r_arrival - r_p
                r = r_p + photoemission_boundary_distance * dvec / norm2(dvec)

                call iu_interpolate_at(pdsim_ug, r, n_vars, &
                     i_vars, vars, i_cell)

                do n = 1, n_photoemission
                   call add_new_avalanche(rng, time, r, vars, avalanches, pq)
                   ! Exit when the inception threshold has been reached
                   if (pq%n_stored == n_stored_limit) exit
                end do
             end if
          end if

          ! Exit when the inception threshold has been reached
          if (pq%n_stored == n_stored_limit) exit
       end do
    else
       ! No photon tracing
       do k = 1, n_photons
          i_cell = 0
          r = pdsim_convert_r(absorption_locations(:, k))
          call iu_get_cell(pdsim_ug, r, i_cell)

          if (i_cell > 0) then
             if (pdsim_ug%icell_data(i_cell, pdsim_icdata_material) == &
                  pdsim_gas_material_value) then
                ! Photoelectron can contribute to new avalanche
                call iu_interpolate_at(pdsim_ug, r, n_vars, &
                     i_vars, vars, i_cell)

                call add_new_avalanche(rng, time, r, vars, avalanches, pq)

                ! Exit when the inception threshold has been reached
                if (pq%n_stored == n_stored_limit) exit
             end if
          end if
       end do
    end if

  end subroutine photoi_SEE

  !> Sample positive ion secondary electron emission from an avalanche and
  !> store the resulting new avalanches
  subroutine ion_SEE(av, time, inception_count, rng, avalanches, pq)
    type(avalanche_t), intent(in)    :: av
    real(dp), intent(in)             :: time
    integer, intent(in)              :: inception_count
    type(rng_t), intent(inout)       :: rng
    type(avalanche_t), intent(inout) :: avalanches(:)
    type(pqr_t), intent(inout)       :: pq

    real(dp) :: mean, r(3), vars(n_vars)
    integer  :: n_secondary_electrons, i_cell, k

    ! Sample secondary emission due to ions
    mean = av%ion_gamma * av%num_ionizations
    n_secondary_electrons = rng%poisson(mean)

    if (n_secondary_electrons > 0) then
       ! Get parameters for avalanches starting at boundary
       i_cell = 0
       call iu_interpolate_at(pdsim_ug, av%r_ion_arrival, &
            n_vars, i_vars, vars, i_cell)

       do k = 1, n_secondary_electrons
          call add_new_avalanche(rng, time + av%ion_travel_time, r, &
               vars, avalanches, pq)

          ! Exit when the inception threshold has been reached
          if (pq%n_stored == inception_count) exit
       end do
    end if
  end subroutine ion_SEE

  !> Add a new avalanche to the list of avalanches, if it has a nonzero size
  subroutine add_new_avalanche(rng, time, r, vars, avalanches, pq, &
       u01_created, u01_size)
    use m_random
    use iso_fortran_env, only: int64
    type(rng_t), intent(inout)       :: rng
    real(dp), intent(in)             :: time
    real(dp), intent(in)             :: r(3)
    real(dp), intent(in)             :: vars(n_vars)
    type(avalanche_t), intent(inout) :: avalanches(:)
    type(pqr_t), intent(inout)       :: pq
    !> Uniform [0, 1) random number, determines if an avalanche is produced
    real(dp), intent(in), optional   :: u01_created
    !> Uniform [0, 1) random number used to sample avalanche size
    real(dp), intent(in), optional   :: u01_size

    integer  :: ix
    real(dp) :: pgeom, unif_01
    real(dp) :: p_m1, k_star, travel_time, r_arrival(3)
    real(dp) :: r_ion_arrival(3), ion_gamma, ion_time
    type(avalanche_t) :: av

    real(dp), parameter :: eps = 1e-12_dp

    ! Unpack vars(:)
    p_m1          = vars(1)
    k_star        = vars(2)
    travel_time   = vars(3)
    r_arrival     = vars(4:6)
    r_ion_arrival = vars(7:9)
    ion_gamma     = vars(10)
    ion_time      = vars(11)

    if (present(u01_created)) then
       unif_01 = u01_created
    else
       unif_01 = rng%unif_01()
    end if

    ! Check if first electron produces additional ionization
    if (unif_01 > p_m1) then

       ! Probability of geometric distribution
       pgeom = (1 - p_m1) / (exp(k_star) - 1)

       if (present(u01_size)) then
          unif_01 = u01_size
       else
          unif_01 = rng%unif_01()
       end if

       ! Sample avalanche size from geometric distribution + 1. Take care of
       ! cases when pgeom >= 1.0 (due to numerical errors) or pgeom ~ 0
       if (pgeom < 1 - eps) then
          av%num_ionizations = geometric(unif_01, pgeom)
       else
          ! exp(k_star) - 1 is small, so avalanche should have zero size
          return
       end if

       ! Add one if it does not cause overflow
       if (av%num_ionizations < huge(1_int64)) then
          av%num_ionizations = av%num_ionizations + 1
       end if

       av%t_source = time
       av%r_source = r
       av%r_arrival = r_arrival
       av%r_ion_arrival = r_ion_arrival
       av%ion_gamma = ion_gamma
       av%ion_travel_time = ion_time
       av%t_arrival = time + travel_time

       ! Get index for avalanche
       call pqr_push_aix(pq, ix, av%t_arrival)

       if (ix > size(avalanches)) then
          print *, "Avalanche index: ", ix
          print *, "Simulation time: ", time
          print *, "Number of active avalanches: ", pq%n_stored
          error stop "Not enough storage for avalanches"
       end if

       ! Store avalanche
       avalanches(ix) = av
    end if

  end subroutine add_new_avalanche

  !> Compute log(1+x) with good accuracy, see "What Every Computer Scientist
  !> Should Know About Floating-Point Arithmetic"
  real(dp) function log1p(x)
    real(dp), intent(in) :: x

    if (1.0_dp + abs(x) > 1.0_dp) then
       log1p = log(1.0_dp + x) * x / ((1.0_dp + x) - 1.0_dp)
    else
       log1p = x
    endif
  end function log1p

  !> Return exponential random variate with a rate of one
  real(dp) function exponential_standard(unif_01)
    real(dp), intent(in) :: unif_01

    if (unif_01 < 0.5_dp) then
       exponential_standard = -log1p(-unif_01)
    else
       exponential_standard = -log(1 - unif_01)
    end if
  end function exponential_standard

  !> Sample from geometric distribution with Pr(X = k) = (1 - p)^(k-1) * p
  integer(int64) function geometric(unif_01, p)
    real(dp), intent(in) :: unif_01
    real(dp), intent(in) :: p
    real(dp)             :: tmp
    real(dp), parameter  :: threshold = real(huge(1_int64) - 1, dp)

    ! Perform inversion sampling X = ceiling(log(U)/log(1-p))
    tmp = -exponential_standard(unif_01) / log1p(-p)

    ! Avoid overflow
    if (tmp < threshold) then
       geometric = ceiling(tmp, int64)
    else
       geometric = huge(1_int64)
    end if

  end function geometric

end module m_avalanche
