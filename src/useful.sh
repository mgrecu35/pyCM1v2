cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding    -DNETCDF module_cm1_vars.F > module_cm1_vars.f90

cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding    -DNETCDF restart.F > restart.f90

cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding    -DNETCDF restart_write.F > restart_write.f90

cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding    -DNETCDF restart_read.F > restart_read.f90

#gfortran -ffree-form -ffree-line-length-none -fPIC -O2 -finline-functions -fopenmp -I/Users/mgrecu/miniconda3//include -c restart.f90

#gfortran -ffree-form -ffree-line-length-none -fPIC -O2 -finline-functions -fopenmp -I/Users/mgrecu/miniconda3//include -c restart_write.f90

gfortran -ffree-form -ffree-line-length-none -fPIC -O2 -finline-functions -fopenmp -I/Users/mgrecu/miniconda3//include -c restart_read.f90
