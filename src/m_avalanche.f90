!> Module for simulation of avalanches
module m_avalanche
  use iso_fortran_env, only: int64
  use m_config
  use m_interp_unstructured
  use m_pdsim
  use m_photoi
  use m_pq

  implicit none
  private

  type avalanche_t
     !> Start time
     real(dp)       :: t_source
     !> Arrival time
     real(dp)       :: t_arrival
     !> Arrival location of avalanche
     real(dp)       :: r_arrival(3)
     !> Arrival location of positive ions
     real(dp)       :: r_ion_arrival(3)
     !> Start location
     real(dp)       :: r_source(3)
     !> Final size of the avalanche
     integer(int64) :: avalanche_size
  end type avalanche_t

  ! Public methods
  public :: avalanche_create_config
  public :: avalanche_simulate

contains

  subroutine avalanche_create_config(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "avalanche%max_total_number", 1000*1000, &
         "Maximum total number of avalanches per run")
    call CFG_add(cfg, "avalanche%n_runs", 10, &
         "Number of runs per initial avalanche location")
    call CFG_add(cfg, "avalanche%max_photons", 10*1000, &
         "Maximum number of photons produced by a single avalanche")
    call CFG_add(cfg, "avalanche%inception_threshold", 1000, &
         "Discharge inception when there are this many future avalanches")

  end subroutine avalanche_create_config

  !> Simulate avalanches as discrete events
  subroutine avalanche_simulate(cfg)
    use m_random
    type(cfg_t), intent(inout) :: cfg

    integer                        :: n, ix, k, i_run, i_avalanche
    integer                        :: max_avalanches, n_runs
    integer                        :: max_photons
    integer                        :: inception_threshold
    integer                        :: n_photons, i_cell, n_secondary_electrons
    integer, parameter             :: n_vars = 9
    integer                        :: i_vars(n_vars)
    real(dp)                       :: r(3), w, k_integral, travel_time
    real(dp)                       :: time, r_arrival(3), r_ion_arrival(3)
    real(dp)                       :: vars(n_vars)
    real(dp)                       :: mean, p_avg, volume
    type(avalanche_t), allocatable :: avalanches(:)
    real(dp), allocatable          :: absorption_locations(:, :)
    real(dp), allocatable          :: inception_time(:)
    logical, allocatable           :: inception(:)
    type(pqr_t)                    :: pq
    type(rng_t)                    :: rng

    i_vars = [i_w, i_k_integral, i_avalanche_time, i_x1, i_x2, i_x3, &
         i_ion_x1, i_ion_x2, i_ion_x3]

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
          travel_time = pdsim_ug%point_data(n, i_avalanche_time)
          r_arrival = pdsim_ug%point_data(n, [i_x1, i_x2, i_x3])
          r_ion_arrival = pdsim_ug%point_data(n, &
               [i_ion_x1, i_ion_x2, i_ion_x3])

          call add_new_avalanche(rng, time, r, w, k_integral, travel_time, &
               r_arrival, r_ion_arrival, i_avalanche, avalanches, pq)

          time_loop: do while (pq%n_stored > 0)
             ! Get the next avalanche
             call pqr_pop(pq, ix, time)

             associate (av => avalanches(ix))
               if (photoi_enabled) then
                  ! Photons are assumed to all originate from the end position
                  ! of the avalanche
                  call photoi_sample_photons(rng, av%r_arrival, &
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
                        r_ion_arrival = vars(7:9)

                        call add_new_avalanche(rng, time, r, w, k_integral, &
                             travel_time, r_arrival, r_ion_arrival, &
                             i_avalanche, avalanches, pq)

                        ! Exit when the inception threshold has been reached
                        if (pq%n_stored == inception_threshold) exit time_loop
                     end if
                  end do
               end if

               if (pdsim_ion_gamma_boundary > 0.0_dp) then
                  ! Sample secondary emission due to ions
                  mean = pdsim_ion_gamma_boundary * av%avalanche_size
                  n_secondary_electrons = rng%poisson(mean)

                  if (n_secondary_electrons > 0) then
                     ! Get parameters for avalanches starting at boundary
                     i_cell = 0
                     call iu_interpolate_at(pdsim_ug, av%r_ion_arrival, &
                          n_vars, i_vars, vars, i_cell)

                     w = vars(1)
                     k_integral = vars(2)
                     travel_time = vars(3)
                     r_arrival = vars(4:6)
                     r_ion_arrival = vars(7:9)

                     do k = 1, n_secondary_electrons
                        call add_new_avalanche(rng, time, r, w, k_integral, &
                             travel_time, r_arrival, r_ion_arrival, &
                             i_avalanche, avalanches, pq)

                        ! Exit when the inception threshold has been reached
                        if (pq%n_stored == inception_threshold) exit time_loop
                     end do
                  end if
               end if
             end associate

          end do time_loop

          inception(i_run) = (pq%n_stored >= inception_threshold)
          inception_time(i_run) = time
       end do

       pdsim_ug%point_data(n, i_inception_time) = &
            sum(inception_time, mask=inception)/max(1, count(inception))
       pdsim_ug%point_data(n, i_inception_prob) = count(inception) / real(n_runs, dp)
    end do

    call pdsim_pointdata_average(pdsim_ug, i_inception_prob, &
         pdsim_axisymmetric, p_avg, volume)

    write(*, "(A,E11.3)") " Average inception probability: ", p_avg
    write(*, "(A,E12.4)") " Total volume of gas: ", volume

  end subroutine avalanche_simulate

  !> Add a new avalanche to the list of avalanches, if it has a nonzero size
  subroutine add_new_avalanche(rng, time, r, w, k_integral, &
       dt, r_arrival, r_ion_arrival, ix, avalanches, pq)
    use m_random
    use iso_fortran_env, only: int64
    type(rng_t), intent(inout)       :: rng
    real(dp), intent(in)             :: time
    real(dp), intent(in)             :: r(3)
    real(dp), intent(in)             :: w
    real(dp), intent(in)             :: k_integral
    real(dp), intent(in)             :: dt
    real(dp), intent(in)             :: r_arrival(3)
    real(dp), intent(in)             :: r_ion_arrival(3)
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
         av%r_ion_arrival = r_ion_arrival
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

end module m_avalanche
