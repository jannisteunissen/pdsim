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
     !> Number of ionizations, including the one that created the first
     !> electron, so this number starts at 1.
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
    call CFG_add(cfg, "avalanche%max_photons", 100*1000, &
         "Maximum number of photons produced by a single avalanche")
    call CFG_add(cfg, "avalanche%inception_count", 1000, &
         "Inception occurs when there are this many future avalanches")
    call CFG_add(cfg, "avalanche%inception_size", 1e8_dp, &
         "Inception occurs when an avalanche of at least this size occurs")
    call CFG_add(cfg, "avalanche%trace_photons", .false., &
         "Whether to trace photon paths")
    call CFG_add(cfg, "avalanche%photoemission_boundary_distance", 1e-9_dp, &
         "Produce photoemission this distance away from the boundary (m)")

  end subroutine avalanche_create_config

  !> Simulate avalanches as discrete events
  subroutine avalanche_simulate(cfg)
    use omp_lib
    type(cfg_t), intent(inout) :: cfg

    integer                        :: n
    integer                        :: n_runs, thread_id
    integer                        :: max_photons
    integer                        :: inception_count
    real(dp)                       :: inception_size
    real(dp)                       :: p_avg, volume
    logical                        :: trace_photons
    type(avalanche_t), allocatable :: avalanches(:)
    real(dp), allocatable          :: inception_time(:)
    logical, allocatable           :: inception(:)
    type(pqr_t)                    :: pq
    type(rng_t)                    :: rng
    type(prng_t)                   :: prng

    call iu_add_point_data(pdsim_ug, "inception_time", i_inception_time)
    call iu_add_point_data(pdsim_ug, "inception_prob", i_inception_prob)

    call CFG_get(cfg, "avalanche%max_photons", max_photons)
    call CFG_get(cfg, "avalanche%n_runs", n_runs)
    call CFG_get(cfg, "avalanche%inception_count", inception_count)
    call CFG_get(cfg, "avalanche%inception_size", inception_size)
    call CFG_get(cfg, "avalanche%trace_photons", trace_photons)
    call CFG_get(cfg, "avalanche%photoemission_boundary_distance", &
         photoemission_boundary_distance)

    allocate(avalanches(inception_count))
    allocate(inception(n_runs))
    allocate(inception_time(n_runs))

    ! Create priority queue that will store upcoming avalanches, sorted by
    ! their arrival time
    call pqr_create(pq, inception_count)

    call rng%set_random_seed()
    call prng%init_parallel(omp_get_max_threads(), rng)

    !$omp parallel private(n, thread_id, pq, avalanches, &
    !$omp &inception, inception_time)

    ! Get a random number generator for each thread
    thread_id = omp_get_thread_num() + 1
    rng = prng%rngs(thread_id)

    !$omp do schedule(dynamic)
    do n = 1, pdsim_ug%n_points
       if (modulo(n, pdsim_ug%n_points/100) == 0) then
          write(*, "(F6.1,A)") (n*1e2_dp)/pdsim_ug%n_points, "%"
       end if

       call run_avalanche(n, n_runs, inception_count, inception_size, &
            max_photons, trace_photons, rng, pq, avalanches, &
            inception, inception_time)

       pdsim_ug%point_data(n, i_inception_time) = &
            sum(inception_time, mask=inception)/max(1, count(inception))
       pdsim_ug%point_data(n, i_inception_prob) = count(inception) / real(n_runs, dp)
    end do
    !$omp end do
    !$omp end parallel

    call pdsim_pointdata_average(pdsim_ug, i_inception_prob, &
         pdsim_axisymmetric, p_avg, volume)

    write(*, "(A,E11.3)") " Average inception probability: ", p_avg
    write(*, "(A,E12.4)") " Total volume of gas: ", volume

  end subroutine avalanche_simulate

  !> Run n_runs avalanches starting at a point
  subroutine run_avalanche(ip, n_runs, inception_count, inception_size, &
       max_photons, trace_photons, rng, pq, avalanches, &
       inception, inception_time)
    integer, intent(in)              :: ip !< Point index
    integer, intent(in)              :: n_runs
    integer, intent(in)              :: inception_count
    real(dp), intent(in)             :: inception_size
    integer, intent(in)              :: max_photons
    logical, intent(in)              :: trace_photons
    type(rng_t), intent(inout)       :: rng
    type(pqr_t), intent(inout)       :: pq
    type(avalanche_t), intent(inout) :: avalanches(inception_count)
    logical, intent(out)             :: inception(n_runs)
    real(dp), intent(out)            :: inception_time(n_runs)

    integer           :: ix, i_run
    real(dp)          :: time, r(3), vars(n_vars)
    type(avalanche_t) :: av

    ! Set module-level variable
    i_vars = [i_p_m1, i_kstar, i_avalanche_time, i_x1, i_x2, i_x3, &
         i_ion_x1, i_ion_x2, i_ion_x3, i_ion_gamma, i_ion_time]

    do i_run = 1, n_runs
       time = 0.0_dp
       call pqr_reset(pq)

       ! Parameters for the initial avalanche
       r = pdsim_ug%points(:, ip)
       vars = pdsim_ug%point_data(ip, i_vars)
       call add_new_avalanche(rng, time, r, vars, avalanches, pq)

       inception(i_run) = .false.
       inception_time(i_run) = 0.0_dp

       do while (pq%n_stored > 0)
          ! Get the next avalanche
          call pqr_pop_aix(pq, ix, time)

          ! Store copy of avalanche, as the index might be overwritten
          av = avalanches(ix)

          ! Check if avalanche size exceeds threshold
          if (avalanches(ix)%num_ionizations > inception_size) then
             inception(i_run) = .true.
             exit
          end if

          if (photoi_enabled) then
             call photoi_SEE(av, time, max_photons, &
                  trace_photons, inception_count, rng, avalanches, pq)

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
          inception_time(i_run) = time
       end if
    end do

  end subroutine run_avalanche

  !> Sample photoionization secondary electron emission from an avalanche and
  !> store the resulting new avalanches
  subroutine photoi_SEE(av, time, max_photons, trace_photons, &
       inception_count, rng, avalanches, pq)
    type(avalanche_t), intent(in)    :: av
    real(dp), intent(in)             :: time
    integer, intent(in)              :: max_photons
    logical, intent(in)              :: trace_photons
    integer, intent(in)              :: inception_count
    type(rng_t), intent(inout)       :: rng
    type(avalanche_t), intent(inout) :: avalanches(:)
    type(pqr_t), intent(inout)       :: pq

    real(dp) :: absorption_locations(3, max_photons)
    integer  :: k, n_photons, i_cell, i_start, status, n_steps, material
    real(dp) :: r(3), r_p(3), vars(n_vars), dvec(3), gamma

    ! Subtract one since initial ionization does not produce photons
    call photoi_sample_photons(rng, av%r_arrival, &
         av%num_ionizations - 1.0_dp, max_photons, &
         n_photons, absorption_locations)

    if (n_photons == 0) return

    if (trace_photons) then
       ! Get start location
       i_start = 0
       call iu_get_cell(pdsim_ug, av%r_arrival, i_start)

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

             if (rng%unif_01() < gamma) then
                ! Create electron due slightly away from boundary
                dvec = av%r_arrival - r_p
                r = r_p + photoemission_boundary_distance * dvec / norm2(dvec)

                call iu_interpolate_at(pdsim_ug, r, n_vars, &
                     i_vars, vars, i_cell)

                call add_new_avalanche(rng, time, r, vars, avalanches, pq)

                ! Exit when the inception threshold has been reached
                if (pq%n_stored == inception_count) exit
             end if
          end if

          ! Exit when the inception threshold has been reached
          if (pq%n_stored == inception_count) exit
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
                if (pq%n_stored == inception_count) exit
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
  subroutine add_new_avalanche(rng, time, r, vars, avalanches, pq)
    use m_random
    use iso_fortran_env, only: int64
    type(rng_t), intent(inout)       :: rng
    real(dp), intent(in)             :: time
    real(dp), intent(in)             :: r(3)
    real(dp), intent(in)             :: vars(n_vars)
    type(avalanche_t), intent(inout) :: avalanches(:)
    type(pqr_t), intent(inout)       :: pq

    integer  :: ix
    real(dp) :: pgeom, tmp
    real(dp) :: p_m1, k_star, travel_time, r_arrival(3)
    real(dp) :: r_ion_arrival(3), ion_gamma, ion_time
    type(avalanche_t) :: av

    ! Unpack vars(:)
    p_m1          = vars(1)
    k_star        = vars(2)
    travel_time   = vars(3)
    r_arrival     = vars(4:6)
    r_ion_arrival = vars(7:9)
    ion_gamma     = vars(10)
    ion_time      = vars(11)

    ! Check if first electron produces additional ionization
    if (rng%unif_01() > p_m1) then

       ! Probability of geometric distribution
       pgeom = (1 - p_m1) / (exp(k_star) - 1)

       ! Sample avalanche size from geometric distribution. Take care of cases
       ! when pgeom >= 1.0 due to numerical errors
       if (pgeom < 1) then
          tmp = 1 + log(1 - rng%unif_01()) / log(1 - pgeom)
       else
          ! exp(k_star) - 1 is small, so avalanche should have zero size
          return
       end if

       av%num_ionizations = ceiling(tmp, int64)
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

end module m_avalanche
