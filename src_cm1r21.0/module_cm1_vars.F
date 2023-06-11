module cm1vars
      use input
      use constants
      use param_module
      integer :: nstep,nstep0
      integer :: srec,sirec,urec,vrec,wrec,nrec,mrec,prec,trecs,trecw,arecs,arecw
      integer :: nrst,nstatout,nwrite,nwritet,nwritea,nwriteh
      integer :: rbufsz,num_soil_layers,ndt
      real :: dt,dtlast
      double precision :: mtime,stattim,taptim,rsttim,radtim,prcltim,diagtim,azimavgtim,hifrqtim
      double precision :: adt,adtlast,acfl,dbldt
      double precision :: mass1,mass2
      double precision :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      logical :: dosfcflx
      logical, dimension(maxq) :: cloudvar,rhovar
      real, dimension(maxq) :: qmag
      character(len=3), dimension(maxq) :: qname
      character(len=20), dimension(maxq) :: qunit
      character(len=6), dimension(maxq) :: budname
      character(len=60), dimension(maxvars) :: desc_output
      character(len=40), dimension(maxvars) :: name_output,unit_output
      character(len=1),  dimension(maxvars) :: grid_output
      logical, dimension(maxvars) :: cmpr_output
      character(len=40), dimension(maxvars) :: name_stat,desc_stat,unit_stat
      character(len=40), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      double precision, dimension(:), allocatable :: bud,bud2
      double precision, dimension(:), allocatable :: qbudget
      double precision, dimension(:), allocatable :: asq,bsq
      real, dimension(:), allocatable :: xh,rxh,arh1,arh2,uh,ruh
      real, dimension(:), allocatable :: xf,rxf,arf1,arf2,uf,ruf
      real, dimension(:), allocatable :: yh,vh,rvh
      real, dimension(:), allocatable :: yf,vf,rvf
      real, dimension(:), allocatable :: xfref,xhref,yfref,yhref
      double precision, dimension(:), allocatable :: dumk1,dumk2,dumk3,dumk4
      double precision, dimension(:,:), allocatable :: dum2d1,dum2d2,dum2d3,dum2d4,dum2d5
      real, dimension(:), allocatable :: rds,sigma,rdsf,sigmaf
      real, dimension(:), allocatable :: wprof,ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,  &
                                         uavg,vavg,savg,thavg,pavg,ulspg,vlspg
      real, dimension(:,:), allocatable :: qavg
      double precision, dimension(:,:), allocatable :: cavg
      real, dimension(:,:,:), allocatable :: tauh,taus,zh,mh,rmh,c1,c2
      real, dimension(:,:,:), allocatable :: tauf,zf,mf,rmf
      real, dimension(:), allocatable :: rstat
      real, dimension(:,:), allocatable :: rho0s,pi0s,prs0s,rth0s
      real, dimension(:,:,:), allocatable :: pi0,rho0,prs0,thv0,th0,rth0,qv0
      real, dimension(:,:,:), allocatable :: qc0,qi0,rr0,rf0,rrf0,u0,v0,thrd
      real, dimension(:,:,:), allocatable :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9
      real, dimension(:,:), allocatable :: zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, dimension(:,:,:), allocatable :: gx,gxu,gy,gyv
      real, dimension(:,:,:), allocatable :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(:,:), allocatable :: tsk,znt,rznt,zntmp,ust,stau,tst,qst,z0t,z0q,thflux,qvflux,  &
                                           cd,ch,cq,u1,v1,s1,t1,xland,psfc,tlh,f2d,psmth,prate,ustt,ut,vt,st,cm0
      real, dimension(:,:), allocatable :: radbcw,radbce
      real, dimension(:,:), allocatable :: radbcs,radbcn
      real, dimension(:,:), allocatable :: dtu,dtu0,dtv,dtv0
      real, dimension(:,:,:), allocatable :: divx,rho,rr,rf,prs
      real, dimension(:,:,:), allocatable :: t11,t12,t13,t22,t23,t33
      real, dimension(:,:,:), allocatable :: m11,m12,m13,m22,m23,m33
      real, dimension(:,:,:), allocatable :: rru,ua,u3d,uten,uten1
      real, dimension(:,:,:), allocatable :: rrv,va,v3d,vten,vten1
      real, dimension(:,:,:), allocatable :: rrw,wa,w3d,wten,wten1
      real, dimension(:,:,:), allocatable :: ppi,pp3d,ppten,sten,sadv,ppx,phi1,phi2
      real, dimension(:,:,:), allocatable :: tha,th3d,thten,thten1,thterm
      real, dimension(:,:,:), allocatable :: qpten,qtten,qvten,qcten
      real, dimension(:,:,:,:), allocatable :: qa,q3d,qten
      real, dimension(:,:,:), allocatable :: p3a
      real, dimension(:,:,:,:), allocatable :: p3o
      real, dimension(:,:,:), allocatable :: kmh,kmv,khh,khv,cme,csm,ce1,ce2
      real, dimension(:,:,:), allocatable :: tkea,tke3d,tketen
      real, dimension(:,:,:), allocatable :: nm,defv,defh,lenscl,dissten,epst,epsd1,epsd2
      real, dimension(:,:,:), allocatable :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,qrpten,qspten,qgpten
      real, dimension(:,:,:), allocatable :: xkzh,xkzq,xkzm
      real, dimension(:,:,:), allocatable :: tsq,qsq,cov,sh3d,el_pbl,qc_bl,qi_bl,cldfra_bl,  &
                                      qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,   &
                                      edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc, &
                                      sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D
      real, dimension(:,:), allocatable :: vdfg,maxmf
      integer, dimension(:,:), allocatable :: nupdraft,ktop_plume
      real, dimension(:,:,:), allocatable :: swten,lwten,swtenc,lwtenc,cldfra,o30
      real, dimension(:,:), allocatable :: zir,radsw,rnflx,radswnet,radlwin,dsr,olr
      real, dimension(:,:,:), allocatable :: rad2d,effc,effi,effs,effr,effg,effis
      real, dimension(:,:), allocatable :: lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc
      real, dimension(:,:), allocatable :: swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc
      real, dimension(:,:), allocatable :: lwcf,swcf,coszr
      real, dimension(:,:), allocatable :: xice,xsnow,xlat,xlong,coszen,swddir,swddni,swddif,hrang
      integer, dimension(:,:,:), allocatable :: cldfra1_flag
      integer, dimension(:,:), allocatable :: lu_index,kpbl2d
      real, dimension(:,:), allocatable :: u10,v10,s10,hfx,qfx,               &
                                      hpbl,wspd,phim,phih,psim,psih,psiq,gz1oz0,br,brcr, &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,       &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,snowh,qcg,dsxy,wstar,delta,prkpp,fm,fh
      real, dimension(:,:), allocatable :: charn,msang,scurx,scury,zkmax,cd_out,ch_out,wscale,wscaleu
      real, dimension(:,:), allocatable :: mznt,swspd,smois,taux,tauy,hpbl2d,evap2d,heat2d

      ! for myj:
      real, dimension(:,:), allocatable :: mixht,akhs,akms,elflx,ct,snow,sice,thz0,qz0,uz0,vz0,u10e,v10e,th10,q10,tshltr,qshltr,pshltr,z0base,zntmyj
      integer, dimension(:,:), allocatable :: lowlyr,ivgtyp
      real, dimension(:,:,:), allocatable :: tke_myj,el_myj
      real, dimension(:,:,:), allocatable :: tmp_pbl

      real, dimension(:), allocatable :: slab_zs,slab_dzs
      real, dimension(:,:,:), allocatable :: tslb
      real, dimension(:,:), allocatable :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, dimension(:,:,:,:),  allocatable :: pta,pt3d,ptten
      real, dimension(:,:), allocatable :: dat1,dat2
      real, dimension(:,:,:), allocatable :: dat3
      integer, dimension(:), allocatable :: reqt
      integer :: reqc,reqk
      real, dimension(:,:), allocatable :: pdata,ploc
      logical, dimension(:,:,:), allocatable :: flag

!--- arrays for MPI ---
      integer, dimension(:), allocatable :: myi1p,myi2p,myj1p,myj2p
      integer, dimension(:), allocatable :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_p2,reqs_p3,reqs_x,reqs_y,reqs_z,reqs_tk
      integer, dimension(:,:),  allocatable :: reqs_q,reqs_t
      real, dimension(:), allocatable :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, dimension(:,:,:), allocatable :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      real, dimension(:,:), allocatable :: ww1,ww2,we1,we2
      real, dimension(:,:), allocatable :: ws1,ws2,wn1,wn2
      real, dimension(:,:), allocatable :: pw1,pw2,pe1,pe2
      real, dimension(:,:), allocatable :: ps1,ps2,pn1,pn2
      real, dimension(:,:), allocatable :: p2w1,p2w2,p2e1,p2e2
      real, dimension(:,:), allocatable :: p2s1,p2s2,p2n1,p2n2
      real, dimension(:,:), allocatable :: p3w1,p3w2,p3e1,p3e2
      real, dimension(:,:), allocatable :: p3s1,p3s2,p3n1,p3n2
      real, dimension(:,:), allocatable :: vw1,vw2,ve1,ve2
      real, dimension(:,:), allocatable :: vs1,vs2,vn1,vn2
      real, dimension(:,:), allocatable :: zw1,zw2,ze1,ze2
      real, dimension(:,:), allocatable :: zs1,zs2,zn1,zn2
      real, dimension(:,:,:), allocatable :: uw31,uw32,ue31,ue32
      real, dimension(:,:,:), allocatable :: us31,us32,un31,un32
      real, dimension(:,:,:), allocatable :: vw31,vw32,ve31,ve32
      real, dimension(:,:,:), allocatable :: vs31,vs32,vn31,vn32
      real, dimension(:,:,:), allocatable :: ww31,ww32,we31,we32
      real, dimension(:,:,:), allocatable :: ws31,ws32,wn31,wn32
      real, dimension(:,:,:), allocatable :: sw31,sw32,se31,se32
      real, dimension(:,:,:), allocatable :: ss31,ss32,sn31,sn32
      real, dimension(:,:,:), allocatable :: rw31,rw32,re31,re32
      real, dimension(:,:,:), allocatable :: rs31,rs32,rn31,rn32
      real, dimension(:,:,:,:), allocatable :: qw31,qw32,qe31,qe32
      real, dimension(:,:,:,:), allocatable :: qs31,qs32,qn31,qn32
      real, dimension(:,:,:), allocatable :: tkw1,tkw2,tke1,tke2
      real, dimension(:,:,:), allocatable :: tks1,tks2,tkn1,tkn2
      real, dimension(:,:,:), allocatable :: kw1,kw2,ke1,ke2
      real, dimension(:,:,:), allocatable :: ks1,ks2,kn1,kn2
      real, dimension(:,:,:,:), allocatable :: tw1,tw2,te1,te2
      real, dimension(:,:,:,:), allocatable :: ts1,ts2,tn1,tn2

      ! arrays for elliptic solver:
      real, dimension(:,:,:),    allocatable :: cfb
      real, dimension(:),        allocatable :: cfa,cfc,d1,d2
      complex, dimension(:,:,:), allocatable :: pdt,lgbth,lgbph
      complex, dimension(:,:),   allocatable :: rhs,trans

      ! diagnostic arrays:
      real, dimension(:,:,:,:), allocatable :: tdiag,qdiag,udiag,vdiag,wdiag,kdiag,pdiag

      ! miscellaneous output:
      real, dimension(:,:,:),   allocatable :: out2d
      real, dimension(:,:,:,:), allocatable :: out3d

      ! for output using cylindrical coords:
      integer, dimension(:,:), allocatable :: cir
      real, dimension(:,:), allocatable :: crr,cangle

      ! t2p arrays:
      real, dimension(:), allocatable :: gamk,gamwall,kmw,ufw,vfw,u1b,v1b,l2p,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4
      real, dimension(:,:,:), allocatable :: u2pt,v2pt,kmwk,ufwk,vfwk

      ! cm1ib arrays:
      logical, dimension(:,:,:), allocatable :: bndy
      integer, dimension(:,:), allocatable :: kbdy
      integer, dimension(:,:,:), allocatable :: hflxw,hflxe,hflxs,hflxn

      logical :: dorestart,dowriteout,dostat,doprclout,dotdwrite,doazimwrite,dohifrqwrite

      ! eddy recycling:
      integer, dimension(:,:), allocatable :: recy_cap,recy_inj
      real, dimension(:,:,:,:), allocatable :: recywe
      real, dimension(:,:,:,:), allocatable :: recysn

      ! for large-scale nudging:
      real, dimension(:,:), allocatable :: lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv

      ! arrays for time-averaging:
      integer, dimension(:), allocatable :: ntavg
      double precision, dimension(:), allocatable :: rtavg
      double precision, dimension(:,:,:,:,:), allocatable :: tavg
      real, dimension(:,:,:,:), allocatable :: timavg
      double precision, dimension(:,:,:,:), allocatable :: sfctavg
      real, dimension(:,:,:), allocatable :: sfctimavg

!-----

      integer count,rate,maxr
      real :: rtime,xtime,time_solve,time_solve0,cfl_target,ks_target,kminit,khinit
      real :: steptime1,steptime2
      integer :: i,j,k,k1,k2,n,nn,fnum,frec
      real :: sum,tem0,tem
      logical :: getsfc,getpbl,update_sfc,startup,restarted,restart_prcl,reset
      logical :: dosolve,dorad,getdbz,getvt,doit,dotbud,doqbud,doubud,dovbud,dowbud,donudge

      integer :: icrs,icenter,jcenter
      real :: xcenter,ycenter
      double precision :: domainlocx,domainlocy

      double precision :: adaptmovetim
      integer :: mvrec,nwritemv

      integer :: ntot
      double precision :: rtot,rdenom
      real(kind=qp) :: qtem
      logical :: printit

#ifdef MPI
      integer :: reqs,rc,ii,jj,id,itmp1,itmp2,jtmp1,jtmp2
      integer, dimension(MPI_STATUS_SIZE) :: status
      real :: mp_total,minval,temx,temy,temni,temnj
      double precision :: tstart,tend
      integer, dimension(:), allocatable :: isum,jsum
#endif
end module cm1vars


module forcing_vars
   real, allocatable :: mufrc(:), mvfrc(:), mqvfrc(:), mthfrc(:)
end module forcing_vars
