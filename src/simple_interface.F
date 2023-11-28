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
     kbo,keo,numqo,thaten_mp_out,th3dten_mp_out,qaten_mp_out,q3dten_mp_out,dt_out,&
     wout,q3dout,th3dout,prsout)
  use cm1vars
  implicit none
  integer :: nsteps, imicro
  integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo
  real, intent(out) :: thaten_mp_out(ibo:ieo,jbo:jeo,kbo:keo)
  real, intent(out) :: th3dten_mp_out(ibo:ieo,jbo:jeo,kbo:keo)
  real, intent(out) :: qaten_mp_out(ibo:ieo,jbo:jeo,kbo:keo,numqo)
  real, intent(out) :: q3dten_mp_out(ibo:ieo,jbo:jeo,kbo:keo,numqo)
  real, intent(out) :: time_out, dt_out
  real, intent(out) :: wout(ibo:ieo,jbo:jeo,kbo:keo+1)
  real, intent(out) :: q3dout(ibo:ieo,jbo:jeo,kbo:keo,numqo)
  real, intent(out) :: th3dout(ibo:ieo,jbo:jeo,kbo:keo)
  real, intent(out) :: prsout(ibo:ieo,jbo:jeo,kbo:keo)
    
  call cm1_timestep(nsteps,time_out,imicro)
  !print*, nsteps,time_out,imicro      
  time_out=mtime
  thaten_mp_out=thaten_mp
  th3dten_mp_out=th3dten_mp
  qaten_mp_out=qaten_mp
  q3dten_mp_out=q3dten_mp
  wout=w3d
  q3dout=q3dold
  th3dout=t3dold
  prsout=prsold    
  dt_out=dbldt    
      
end subroutine pytimestep

subroutine get_init_cond(ibo,ieo,jbo,jeo,&
     kbo,keo,numqo,t3d_out,q3d_out,w3d_out,ppi_out)
  use cm1vars
  implicit none
  integer :: nsteps, imicro
  integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo
  real, intent(out) :: q3d_out(ibo:ieo,jbo:jeo,kbo:keo,numqo)
  real, intent(out) :: t3d_out(ibo:ieo,jbo:jeo,kbo:keo)
  real, intent(out) :: w3d_out(ibo:ieo,jbo:jeo,kbo:keo+1)    
  real, intent(out) :: ppi_out(ibo:ieo,jbo:jeo,kbo:keo)
  t3d_out=(th3d)
  q3d_out=(q3d)
  w3d_out=w3d
  ppi_out=pp3d   
      
end subroutine get_init_cond

subroutine set_init_cond(ibo,ieo,jbo,jeo,&
     kbo,keo,numqo,t3d_out,q3d_out,w3d_out,u3d_out,v3d_out,ppi_out)
  use cm1vars
  implicit none
  integer :: nsteps, imicro
  integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo
  real  :: q3d_out(ibo:ieo,jbo:jeo,kbo:keo,numqo)
  real  :: t3d_out(ibo:ieo,jbo:jeo,kbo:keo)
  real  :: w3d_out(ibo:ieo,jbo:jeo,kbo:keo+1)
  real  :: u3d_out(ibo:ieo+1,jbo:jeo,kbo:keo)
  real  :: v3d_out(ibo:ieo,jbo:jeo+1,kbo:keo)    
  real  :: ppi_out(ibo:ieo,jbo:jeo,kbo:keo)
  th3d=t3d_out
  tha=t3d_out    
  q3d=q3d_out
  qa=q3d_out
  w3d=w3d_out
  wa=w3d_out
  u3d=u3d_out
  ua=u3d_out
  v3d=v3d_out
  va=v3d_out    
  pp3d=ppi_out
  ppi=ppi_out    
end subroutine set_init_cond

subroutine make_q_nonneg()
  use cm1vars
  implicit none

  integer :: l
  do i=ib,ie
     do j=jb,je
        do k=ke,kb
           do l=1,6
              if(q3d(i,j,k,l)<0) q3d(i,j,k,l)=0
              if(qa(i,j,k,l)<0) qa(i,j,k,l)=0
           enddo
         enddo
     enddo
  enddo
      
  
end subroutine make_q_nonneg
            
      
      
subroutine get_mp_aux(ibo,ieo,jbo,jeo,&
     kbo,keo,numqo,prs_mp_out,t3d_mp_out,ppi3d)
  use cm1vars
  implicit none
  integer :: nsteps, imicro
  integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo
  real, intent(out) :: prs_mp_out(ibo:ieo,jbo:jeo,kbo:keo)
  real, intent(out) :: t3d_mp_out(ibo:ieo,jbo:jeo,kbo:keo)
  real, intent(out) :: ppi3d(ibo:ieo,jbo:jeo,kbo:keo)
  t3d_mp_out=(th0+th3d)*(pi0+pp3d)
  ppi3d=(pi0+pp3d)
  prs_mp_out=prs    
      
end subroutine get_mp_aux

subroutine get_t0(ibo,ieo,jbo,jeo,&
     kbo,keo,numqo,t3d_out,prs0_out,qv0_out,pi0_out)
  use cm1vars
  implicit none
  integer :: nsteps, imicro
  integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo
  real, intent(out) :: t3d_out(ibo:ieo,jbo:jeo,kbo:keo),&
      prs0_out(ibo:ieo,jbo:jeo,kbo:keo),qv0_out(ibo:ieo,jbo:jeo,kbo:keo),&
      pi0_out(ibo:ieo,jbo:jeo,kbo:keo)
  t3d_out=(th0)
  pi0_out=(pi0)
  prs0_out=prs0
  qv0_out=qv0    
end subroutine get_t0
      

      
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


subroutine py_mp(ibo,ieo,jbo,jeo,&
      kbo,keo,numqo,nzp,q3dp,t3dp,q3dtenp,t3dtenp,prsp,dtp,dzp)
use cm1vars
use MODULE_MP_MORR_TWO_MOMENT2      
implicit none      
integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo,nzp
real, intent(in) ::q3dp(ibo:ieo,jbo:jeo,kbo:keo,numqo)
real, intent(in) ::t3dp(ibo:ieo,jbo:jeo,kbo:keo),prsp(ibo:ieo,jbo:jeo,kbo:keo)

real, intent(out) ::q3dtenp(ibo:ieo,jbo:jeo,kbo:keo,numqo)
real, intent(out) ::t3dtenp(ibo:ieo,jbo:jeo,kbo:keo)

integer :: iqv,iqc,iqr,iqi,iqs,iqg,ini,ins,inr,ing,inc
integer :: i1,j1,kb1
real    :: dtp, dzq(nzp),dzp      
iqv=1
iqc=2
iqr=3
iqi=4
iqs=5
iqg=6
ini=7
ins=8
inr=9
ing=10
inc=11
kb1=1
q3dtenp=q3dp
t3dtenp=t3dp
dzq=dzp      
do i1=ibo+3,ieo-3
 do j1=jbo+3,jeo-3     
      call SIMPLE_MICRO3_2(q3dtenp(i1,j1,1:nzp,iqc),q3dtenp(i1,j1,1:nzp,iqi),&
      q3dtenp(i1,j1,1:nzp,iqs),q3dtenp(i1,j1,1:nzp,iqr),&
      t3dtenp(i1,j1,1:nzp),&
      q3dtenp(i1,j1,1:nzp,iqv),prsp(i1,j1,1:nzp),dzq(1:nzp),   &
      dtp,                                                   &
      kb1,nzp,q3dtenp(i1,j1,1:nzp,iqg))
  enddo
enddo
q3dtenp=q3dtenp-q3dp
q3dtenp(:,:,:,7:10)=0      
t3dtenp=(t3dtenp-t3dp)/(pi0+pp3d)      
end subroutine py_mp


subroutine micro_init()
      use MODULE_MP_MORR_TWO_MOMENT2
      call MORR_TWO_MOMENT_INIT2
end subroutine micro_init


subroutine py_qvs(ibo,ieo,jbo,jeo,&
      kbo,keo,numqo,nzp,q3dp,t3dp,prsp,dtp,dzp,qvs3d,xxlv3d,xxls3d,cpm3d)
use cm1vars
use MODULE_MP_MORR_TWO_MOMENT2, only:getqvs      
implicit none      
integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo,nzp
real, intent(in) ::q3dp(ibo:ieo,jbo:jeo,kbo:keo,numqo)
real, intent(in) ::t3dp(ibo:ieo,jbo:jeo,kbo:keo),prsp(ibo:ieo,jbo:jeo,kbo:keo)

real, intent(out) ::qvs3d(ibo:ieo,jbo:jeo,kbo:keo)
real, intent(out) :: xxls3d(ibo:ieo,jbo:jeo,kbo:keo)
real, intent(out) :: xxlv3d(ibo:ieo,jbo:jeo,kbo:keo)
real, intent(out) :: cpm3d(ibo:ieo,jbo:jeo,kbo:keo)
      

integer :: iqv,iqc,iqr,iqi,iqs,iqg,ini,ins,inr,ing,inc
integer :: i1,j1,kb1,kk1
real    :: dtp, dzq(nzp),dzp, xxls(nzp), xxlv(nzp), cpm(nzp)      
iqv=1
iqc=2
iqr=3
iqi=4
iqs=5
iqg=6
ini=7
ins=8
inr=9
ing=10
inc=11
kb1=1
!q3dtenp=q3dp
!      t3dtenp=t3dp
!      (T3D,QV3D,QVS,PRES,KTS,KTE)
dzq=dzp      
do i1=ibo+3,ieo-3
 do j1=jbo+3,jeo-3     
      call getqvs(t3dp(i1,j1,1:nzp),q3dp(i1,j1,1:nzp,iqv),&
      qvs3d(i1,j1,1:nzp),prsp(i1,j1,1:nzp),&
      xxls,xxlv,cpm,kb1,nzp)
      
      xxls3d(i1,j1,1:nzp)=xxls
      xxlv3d(i1,j1,1:nzp)=xxlv
      cpm3d(i1,j1,1:nzp)=cpm
  enddo
enddo
end subroutine py_qvs


subroutine py_autosed(ibo,ieo,jbo,jeo,&
      kbo,keo,numqo,nzp,dt,dzp)
use cm1vars, only: th0, th3d, pi0, pp3d, prs, qa,q3d,rho
use MODULE_MP_MORR_TWO_MOMENT2, only:autoconv_sed      
implicit none
      
integer :: ibo,ieo,jbo,jeo,kbo,keo,numqo,nzp

real :: dzq(1:nzp),rho3d(ibo:ieo,jbo:jeo,kbo:keo)
real :: dt,dzp
integer :: iqv,iqc,iqr,iqi,iqs,iqg,ini,ins,inr,ing,inc
integer :: i1,j1,kb1,kk1

iqv=1
iqc=2
iqr=3
iqi=4
iqs=5
iqg=6
ini=7
ins=8
inr=9
ing=10
inc=11
kb1=1
!print*, rho(1,1,1:nzp)
!stop
dzq=dzp      
do i1=ibo+3,ieo-3
 do j1=jbo+3,jeo-3     
      call autoconv_sed(kb1,nzp,q3d(i1,j1,1:nzp,iqc),q3d(i1,j1,1:nzp,iqr),&
      q3d(i1,j1,1:nzp,iqi),q3d(i1,j1,1:nzp,iqs),rho(i1,j1,1:nzp),dzq,dt)
 enddo
enddo
qa=q3d      
end subroutine py_autosed

subroutine set_prs(ibmo,iemo,jbmo,jemo,kbmo,kemo,prs_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: prs_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
prs=prs_in
end subroutine set_prs

subroutine set_rho(ibmo,iemo,jbmo,jemo,kbmo,kemo,rho_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: rho_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
rho=rho_in
end subroutine set_rho


subroutine get_saturation3d(ibo,ieo,jbo,jeo,kbo,keo,&
      qvs3d,t3d_in,prs_in)
use cm1vars
implicit none      
real :: t3d1,pres1,qvs1,qvi1
integer :: ibo,ieo,jbo,jeo,kbo,keo      
real, intent(in) ::t3d_in(ibo:ieo,jbo:jeo,kbo:keo),&
      prs_in(ibo:ieo,jbo:jeo,kbo:keo)

real, intent(out) ::qvs3d(ibo:ieo,jbo:jeo,kbo:keo)


 do i=1,ieo-3
   do j=1,jeo-3
      do k=1,keo-1
          call saturation_new(t3d_in(i,j,k),prs_in(i,j,k),qvs3d(i,j,k),qvi1)
      enddo
   enddo
 enddo
      
 !     print*, maxval(t3d_in)
 !     print*, maxval(prs_in)
 !     print*, maxval(qvs3d)
end subroutine get_saturation3d
!

      
