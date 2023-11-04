subroutine pycm1_init(ibo,ieo,jbo,jeo,kbo,keo,&
     ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo)
  use cm1vars
  use input
  use forcing_vars    
  implicit none
  integer, intent(out) :: ibo,ieo,jbo,jeo,kbo,keo,&
       ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
 
  call cm1_init()
  ibo=ib
  ieo=ie
  jbo=jb
  jeo=je
  kbo=kb
  keo=ke
  ibmo=ibm
  iemo=iem
  jbmo=jbm
  jemo=jem
  kbmo=kbm
  kemo=kem
  numqo=numq
  allocate(mthfrc(ke-kb-1))
  allocate(mqvfrc(ke-kb-1))
  allocate(mufrc(ke-kb-1))
  allocate(mvfrc(ke-kb-1))
end subroutine pycm1_init

subroutine pytimestep(nsteps,time_out,imicro)
  integer :: nsteps, imicro
  real, intent(out) :: time_out
  call cm1_timestep(nsteps,time_out,imicro)
end subroutine pytimestep
