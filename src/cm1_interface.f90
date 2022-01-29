subroutine pycm1_init(ibo,ieo,jbo,jeo,kbo,keo,&
     ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo)
  use cm1vars
  use input
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
end subroutine pycm1_init

subroutine pytimestep(nsteps,time_out)
  integer :: nsteps
  real, intent(out) :: time_out
  call cm1_timestep(nsteps,time_out)
end subroutine pytimestep


subroutine py_endcm1()
  call end_simulation
end subroutine py_endcm1

subroutine get_tha(ibmo,iemo,jbmo,jemo,kbmo,kemo,th0_out)
  use cm1vars
  implicit none
  integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
  real, intent(out) :: th0_out(ibmo:iemo,jbmo:jemo,kbmo:kemo)
  th0_out=tha
end subroutine get_tha


subroutine set_tha(ibmo,iemo,jbmo,jemo,kbmo,kemo,th0_in)
  use cm1vars
  implicit none
  integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
  real, intent(in) :: th0_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
  tha=th0_in
end subroutine set_tha


subroutine get_th3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,th0_out)
  use cm1vars
  implicit none
  integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
  real, intent(out) :: th0_out(ibmo:iemo,jbmo:jemo,kbmo:kemo)
  th0_out=th3d
end subroutine get_th3d


subroutine set_th3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,th0_in)
  use cm1vars
  implicit none
  integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
  real, intent(in) :: th0_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
  th3d=th0_in
end subroutine set_th3d

!-----------humidity-----------------!

subroutine get_qa(ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo,q0_out)
  use cm1vars
  implicit none
  integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
  real, intent(out) :: q0_out(ibmo:iemo,jbmo:jemo,kbmo:kemo,numqo)
  q0_out=qa
end subroutine get_qa


subroutine set_qa(ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo,q0_in)
  use cm1vars
  implicit none
  integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
  real, intent(in) :: q0_in(ibmo:iemo,jbmo:jemo,kbmo:kemo,numqo)
  qa=q0_in
end subroutine set_qa


subroutine get_q3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo,q0_out)
  use cm1vars
  implicit none
  integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
  real, intent(out) :: q0_out(ibmo:iemo,jbmo:jemo,kbmo:kemo,numqo)
  q0_out=q3d
end subroutine get_q3d


subroutine set_q3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo,q0_in)
  use cm1vars
  implicit none
  integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo,numqo
  real, intent(in) :: q0_in(ibmo:iemo,jbmo:jemo,kbmo:kemo,numqo)
  q3d=q0_in
end subroutine set_q3d

subroutine set_frc(nz1,th_frc,qv_frc)
  use cm1vars
  implicit none
  integer :: nz1
  real, intent(in) :: th_frc(nz1), qv_frc(nz1)
  thfrc(1:nz1)=th_frc
  qvfrc(1:nz1)=qv_frc
end subroutine set_frc
