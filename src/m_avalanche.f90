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
     !> Number of ionizations, including the one that created the first
     !> electron, so this number starts at 1.
     integer(int64) :: num_ionizations
  end type avalanche_t

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
    type(avalanche_t), allocatable :: avalanches(:)
    real(dp), allocatable          :: absorption_locations(:, :)
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

    allocate(avalanches(inception_count))
    allocate(absorption_locations(3, max_photons))
    allocate(inception(n_runs))
    allocate(inception_time(n_runs))

    ! Create priority queue that will store upcoming avalanches, sorted by
    ! their arrival time
    call pqr_create(pq, inception_count)

    call rng%set_random_seed()
    call prng%init_parallel(omp_get_max_threads(), rng)

    !$omp parallel private(n, thread_id, pq, avalanches, absorption_locations, &
    !$omp &inception, inception_time)

    ! Get a random number generator for each thread
    thread_id = omp_get_thread_num() + 1
    rng = prng%rngs(thread_id)

    !$omp do schedule(dynamic)
    do n = 1, pdsim_ug%n_points
       if (modulo(n, pdsim_ug%n_points/100) == 0) then
          write(*, "(F6.1,A)") (n*1e2_dp)/pdsim_ug%n_points, "%"
       end if

       call run_avalanche(n, n_runs, inception_count, inception_size, rng, &
            pq, avalanches, max_photons, absorption_locations, &
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
  subroutine run_avalanche(ip, n_runs, inception_count, inception_size, rng, &
       pq, avalanches, max_photons, absorption_locations, &
       inception, inception_time)
    integer, intent(in)              :: ip !< Point index
    integer, intent(in)              :: n_runs
    integer, intent(in)              :: inception_count
    real(dp), intent(in)             :: inception_size
    type(rng_t), intent(inout)       :: rng
    type(pqr_t), intent(inout)       :: pq
    type(avalanche_t), intent(inout) :: avalanches(inception_count)
    integer, intent(in)              :: max_photons
    real(dp), intent(inout)          :: absorption_locations(3, max_photons)
    logical, intent(out)             :: inception(n_runs)
    real(dp), intent(out)            :: inception_time(n_runs)

    integer, parameter :: n_vars = 10
    integer            :: i_vars(n_vars)
    integer            :: ix, k, i_run, i_cell
    integer            :: n_photons, n_secondary_electrons
    real(dp)           :: r(3), p_m1, k_star, travel_time
    real(dp)           :: time, r_arrival(3), r_ion_arrival(3)
    real(dp)           :: ion_gamma
    real(dp)           :: vars(n_vars), mean

    i_vars = [i_p_m1, i_kstar, i_avalanche_time, i_x1, i_x2, i_x3, &
             i_ion_x1, i_ion_x2, i_ion_x3, i_ion_gamma]

    do i_run = 1, n_runs
       time = 0.0_dp
       call pqr_reset(pq)

       ! Parameters for the initial avalanche
       r = pdsim_ug%points(:, ip)
       p_m1 = pdsim_ug%point_data(ip, i_p_m1)
       k_star = pdsim_ug%point_data(ip, i_kstar)
       travel_time = pdsim_ug%point_data(ip, i_avalanche_time)
       r_arrival = pdsim_ug%point_data(ip, [i_x1, i_x2, i_x3])
       r_ion_arrival = pdsim_ug%point_data(ip, &
            [i_ion_x1, i_ion_x2, i_ion_x3])
       ion_gamma = pdsim_ug%point_data(ip, i_ion_gamma)

       call add_new_avalanche(rng, time, r, p_m1, k_star, travel_time, &
            r_arrival, r_ion_arrival, ion_gamma, avalanches, pq)

       inception(i_run) = .false.

       time_loop: do while (pq%n_stored > 0)
          ! Get the next avalanche
          call pqr_pop_aix(pq, ix, time)

          associate (av => avalanches(ix))

            if (av%num_ionizations > inception_size) then
               inception(i_run) = .true.
               inception_time(i_run) = time
               exit time_loop
            end if

            if (photoi_enabled) then
               ! Subtract one since initial ionization does not produce photons
               call photoi_sample_photons(rng, av%r_arrival, &
                    av%num_ionizations - 1.0_dp, max_photons, &
                    n_photons, absorption_locations)

               do k = 1, n_photons
                  i_cell = 0
                  r = pdsim_convert_r(absorption_locations(:, k))
                  call iu_get_cell(pdsim_ug, r, i_cell)

                  if (i_cell > 0) then
                     ! Photoelectron can contribute to new avalanche
                     call iu_interpolate_at(pdsim_ug, r, n_vars, &
                          i_vars, vars, i_cell)

                     p_m1 = vars(1)
                     k_star = vars(2)
                     travel_time = vars(3)
                     r_arrival = vars(4:6)
                     r_ion_arrival = vars(7:9)
                     ion_gamma = vars(10)

                     call add_new_avalanche(rng, time, r, p_m1, k_star, &
                          travel_time, r_arrival, r_ion_arrival, ion_gamma, &
                          avalanches, pq)

                     ! Exit when the inception threshold has been reached
                     if (pq%n_stored == inception_count) then
                        inception(i_run) = .true.
                        inception_time(i_run) = time
                        exit time_loop
                     end if

                     ! TODO: if no additional ionization was produced, we
                     ! ignore the avalanche, but the probability of secondary
                     ! emission from the single positive ion (from the photon)
                     ! could be included.
                  end if
               end do
            end if

            if (pdsim_ion_see_enabled) then
               ! Sample secondary emission due to ions
               mean = av%ion_gamma * av%num_ionizations
               n_secondary_electrons = rng%poisson(mean)

               if (n_secondary_electrons > 0) then
                  ! Get parameters for avalanches starting at boundary
                  i_cell = 0
                  call iu_interpolate_at(pdsim_ug, av%r_ion_arrival, &
                       n_vars, i_vars, vars, i_cell)

                  p_m1 = vars(1)
                  k_star = vars(2)
                  travel_time = vars(3)
                  r_arrival = vars(4:6)
                  r_ion_arrival = vars(7:9)
                  ion_gamma = vars(10)

                  do k = 1, n_secondary_electrons
                     call add_new_avalanche(rng, time, r, p_m1, k_star, &
                          travel_time, r_arrival, r_ion_arrival, ion_gamma, &
                          avalanches, pq)

                     ! Exit when the inception threshold has been reached
                     if (pq%n_stored == inception_count) then
                        inception(i_run) = .true.
                        inception_time(i_run) = time
                        exit time_loop
                     end if
                  end do
               end if
            end if
          end associate

       end do time_loop
    end do
  end subroutine run_avalanche

  !> Add a new avalanche to the list of avalanches, if it has a nonzero size
  subroutine add_new_avalanche(rng, time, r, p_m1, k_star, &
       dt, r_arrival, r_ion_arrival, ion_gamma, avalanches, pq)
    use m_random
    use iso_fortran_env, only: int64
    type(rng_t), intent(inout)       :: rng
    real(dp), intent(in)             :: time
    real(dp), intent(in)             :: r(3)
    real(dp), intent(in)             :: p_m1
    real(dp), intent(in)             :: k_star
    real(dp), intent(in)             :: dt
    real(dp), intent(in)             :: r_arrival(3)
    real(dp), intent(in)             :: r_ion_arrival(3)
    real(dp), intent(in)             :: ion_gamma
    type(avalanche_t), intent(inout) :: avalanches(:)
    type(pqr_t), intent(inout)       :: pq

    integer  :: ix
    real(dp) :: pgeom, tmp, t_arrival

    ! Check if first electron produces additional ionization
    if (rng%unif_01() > p_m1) then
       ! Get index for avalanche
       t_arrival = time + dt
       call pqr_push_aix(pq, ix, t_arrival)

       if (ix > size(avalanches)) then
          print *, "Avalanche index: ", ix
          print *, "Simulation time: ", time
          print *, "Number of active avalanches: ", pq%n_stored
          error stop "Not enough storage for avalanches"
       end if

       associate (av => avalanches(ix))
         ! Probability of geometric distribution
         pgeom = (1 - p_m1) / (exp(k_star) - p_m1)

         ! Sample avalanche size. Take care of cases when pgeom >= 1.0 due
         ! to numerical errors
         if (pgeom < 1) then
            tmp = log(1 - rng%unif_01()) / log(1 - pgeom)
         else
            tmp = 1.0_dp
         end if

         av%num_ionizations = ceiling(tmp, int64)
         av%t_source = time
         av%r_source = r
         av%r_arrival = r_arrival
         av%r_ion_arrival = r_ion_arrival
         av%ion_gamma = ion_gamma
         av%t_arrival = t_arrival
       end associate
    end if

  end subroutine add_new_avalanche

end module m_avalanche
