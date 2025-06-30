# This file is generated automatically. DO NOT EDIT!

build/pdsim :  \
	build/m_avalanche.o \
	build/m_config.o \
	build/m_integral.o \
	build/m_lookup_table.o \
	build/m_particles.o \
	build/m_pdsim.o \
	build/m_photoi.o \
	build/m_pq.o \
	build/pdsim.o

build/m_avalanche.o :  \
	build/m_config.o \
	build/m_pdsim.o \
	build/m_photoi.o \
	build/m_pq.o

build/m_config.o : 

build/m_integral.o :  \
	build/m_config.o \
	build/m_lookup_table.o \
	build/m_pdsim.o

build/m_lookup_table.o : 

build/m_particles.o :  \
	build/m_config.o \
	build/m_pdsim.o \
	build/m_photoi.o

build/m_pdsim.o :  \
	build/m_config.o \
	build/m_lookup_table.o

build/m_photoi.o :  \
	build/m_config.o

build/m_pq.o : 

build/pdsim.o :  \
	build/m_avalanche.o \
	build/m_config.o \
	build/m_integral.o \
	build/m_particles.o \
	build/m_pdsim.o \
	build/m_photoi.o
