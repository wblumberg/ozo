mod_poisson_interp.o : src/mod_poisson_interp.f90 
mod_wrf_file.o : src/mod_wrf_file.f90 mod_const.o 
mod_const.o : src/mod_const.f90 
ozo.o : src/ozo.f90 mod_time_step_loop.o mod_wrf_file.o mod_const.o 
mod_poisson_DFT.o : src/mod_poisson_DFT.f90 mod_poisson_interp.o mkl_poisson.o mod_poisson_green.o 
mod_time_step_loop.o : src/mod_time_step_loop.f90 mod_common_subrs.o mod_const.o mod_omega.o mod_zo.o mod_wrf_file.o 
mod_zo.o : src/mod_zo.f90 mod_poisson_DFT.o mod_common_subrs.o mod_wrf_file.o mod_const.o 
mod_omega.o : src/mod_omega.f90 mod_common_subrs.o mod_wrf_file.o mod_const.o 
mod_common_subrs.o : src/mod_common_subrs.f90 mod_const.o mod_poisson_DFT.o mod_wrf_file.o 
mod_poisson_green.o : src/mod_poisson_green.f90 
mkl_poisson.o : /home/gblumberg/intel/oneapi/mkl/latest//include/mkl_poisson.f90 mkl_dfti.o 
mkl_dfti.o : /home/gblumberg/intel/oneapi/mkl/latest//include/mkl_dfti.f90 
