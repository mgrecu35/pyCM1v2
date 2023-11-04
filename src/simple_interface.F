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

subroutine pytimestep(nsteps,time_out,imicro,ibo,ieo,jbo,jeo,&
     kbo,keo,numqo,thaten_mp_out,th3dten_mp_out,qaten_mp_out,q3dten_mp_out,dt_out)
  use cm1vars
  implicit none
  integer :: nsteps, imicro
  integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo
  real, intent(out) :: thaten_mp_out(ibo:ieo,jbo:jeo,kbo:keo)
  real, intent(out) :: th3dten_mp_out(ibo:ieo,jbo:jeo,kbo:keo)
  real, intent(out) :: qaten_mp_out(ibo:ieo,jbo:jeo,kbo:keo,numqo)
  real, intent(out) :: q3dten_mp_out(ibo:ieo,jbo:jeo,kbo:keo,numqo)
  real, intent(out) :: time_out, dt_out
  call cm1_timestep(nsteps,time_out,imicro)
  time_out=mtime
  thaten_mp_out=thaten_mp
  th3dten_mp_out=th3dten_mp
  qaten_mp_out=qaten_mp
  q3dten_mp_out=q3dten_mp
  dt_out=dbldt    
      
end subroutine pytimestep


subroutine set_th3d_ten(nxi,nyi,nzi,th0ten,dtp,dt0p)
  use cm1vars
  implicit none
  integer :: nxi,nyi,nzi
  real, intent(in) :: th0ten(nxi,nyi,nzi)
  real:: dtp, dt0p
  th3d(1:nxi,1:nyi,1:nzi)=th3d(1:nxi,1:nyi,1:nzi)+th0ten*dtp/dt0p
  tha(1:nxi,1:nyi,1:nzi)=tha(1:nxi,1:nyi,1:nzi)+th0ten*dtp/dt0p
end subroutine set_th3d_ten

!-----------humidity-----------------!

subroutine set_q3d_ten(nxi,nyi,nzi,numq0,qtenp,dtp,dt0p)
  use cm1vars
  implicit none
  integer :: nxi,nyi,nzi,numq0
  integer :: ip,jp,kp,lp
  real, intent(in) :: qtenp(nxi,nyi,nzi,numq0)
  real:: dtp, dt0p
  q3d(1:nxi,1:nyi,1:nzi,:)=q3d(1:nxi,1:nyi,1:nzi,:)+qtenp*dtp/dt0p
  qa(1:nxi,1:nyi,1:nzi,:)=qa(1:nxi,1:nyi,1:nzi,:)+qtenp*dtp/dt0p
  do ip=1,nxi
     do jp=1,nyi
        do kp=1,nzi
           do lp=1,numq0
              if(q3d(ip,jp,kp,lp)<0) q3d(ip,jp,kp,lp)=0
              if(qa(ip,jp,kp,lp)<0) qa(ip,jp,kp,lp)=0
           enddo
        enddo
    enddo
  enddo
end subroutine set_q3d_ten
