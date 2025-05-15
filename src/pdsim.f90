program simulate_inception
  use m_config
  use m_pdsim
  use m_cross_sec
  use m_particle_core
  use m_units_constants
  use m_gas
  use m_photoi
  use m_lookup_table

  implicit none

  type(CFG_t)   :: cfg
  type(pdsim_t) :: pd

  integer  :: iteration
  real(dp) :: dt, time
  logical  :: output_particles_now
  integer  :: max_iterations = huge(1)
  integer  :: n_output, run_idx

  call pdsim_create_config(cfg)

  call CFG_update_from_arguments(cfg)
  call CFG_check(cfg)

  call pdsim_initialize(cfg, pd)

  ! Repeat the complete simulation "n_runs" times
  do run_idx = 1, pd%n_runs

     call pd%pc%remove_particles()

     call pdsim_create_initial_electron(pd)

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
           call pdsim_write_particles(pd, n_output, run_idx)
        end if

        call pd%pc%advance(dt)
        call pd%pc%after_mover(dt)
        call pd%pc%clean_up()
        time = time + dt

        call pdsim_handle_events(pd%pc)

        if (pd%pc%n_part == 0) exit

        if (pd%pc%get_num_real_part() > pd%num_electrons_inception) exit
     end do

  end do

end program simulate_inception
