program pdsim
  use m_config
  use m_interp_unstructured

  use m_pdsim
  use m_photoi
  use m_integral
  use m_avalanche
  use m_particles

  implicit none

  type(CFG_t) :: cfg

  character(len=40) :: what_to_simulate

  call CFG_add(cfg, "simulate", undefined_str, &
         "What to simulate (integral, avalanches, particles)", &
         required=.true.)

  call pdsim_create_config(cfg)
  call photoi_create_cfg(cfg)
  call integral_create_config(cfg)
  call avalanche_create_config(cfg)
  call particles_create_config(cfg)

  call CFG_update_from_arguments(cfg)
  call CFG_check(cfg)

  call pdsim_initialize(cfg)
  call photoi_initialize(cfg)
  call CFG_get(cfg, "simulate", what_to_simulate)

  ! Always compute ionization integrals
  if (pdsim_verbosity > 0) print *, "Computing ionization integrals"
  call integral_compute(cfg)

  select case (what_to_simulate)
     case ("integral")
        ! Already done above
        continue
     case ("avalanches")
        if (pdsim_verbosity > 0) print *, "Starting avalanche simulation"
        call avalanche_simulate(cfg)
     case ("particles")
        if (pdsim_verbosity > 0) print *, "Starting particle simulation"
        call particles_simulate(cfg)
     case default
        error stop "simulate option can be integral, avalanches, particles"
  end select

  if (pdsim_output_level > 1) then
     call iu_write_vtk(pdsim_ug, trim(pdsim_output_name) // ".vtu")
     if (pdsim_verbosity > 0) &
          print *, "Wrote ", trim(pdsim_output_name) // ".vtu"
  end if

end program pdsim
