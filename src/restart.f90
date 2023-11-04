  MODULE restart_module

  implicit none

  private
  public :: write_restart,read_restart

  CONTAINS

      subroutine write_restart(nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,      &
                               trecs,trecw,arecs,arecw,                             &
                               nwrite,nwritet,nwritea,nwriteh,nrst,nstatout,        &
                               num_soil_layers,nrad2d,                              &
                               avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp, &
                               dt,dtlast,mtime,ndt,adt,acfl,dbldt,mass1,            &
                               stattim,taptim,rsttim,radtim,prcltim,                &
                               qbudget,asq,bsq,qname,                               &
                               xfref,yfref,zh,zf,sigma,sigmaf,zs,                   &
                               th0,prs0,pi0,rho0,qv0,u0,v0,                         &
                               rain,sws,svs,sps,srs,sgs,sus,shs,                    &
                               tsk,znt,ust,cd,ch,cq,u1,v1,s1,t1,thflux,qvflux,      &
                               prate,ustt,ut,vt,st,                                 &
                               radbcw,radbce,radbcs,radbcn,                         &
                               rho,prs,ua,va,wa,ppi,tha,qa,tkea,                    &
                               swten,lwten,radsw,rnflx,radswnet,radlwin,rad2d,      &
                               effc,effi,effs,effr,effg,effis,                      &
                               lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,xland,      &
                               hpbl,wspd,psim,psih,gz1oz0,br,                       &
                               CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                       &
                               MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                    &
                               CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,       &
                               gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta,tslb,    &
                               tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,               &
                               qpten,qtten,qvten,qcten,pta,pdata,ploc,ppx,          &
                               tdiag,qdiag,phi1,phi2,                               &
                   tsq,qsq,cov,sh3d,el_pbl,qc_bl,qi_bl,cldfra_bl,            &
                   qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,            &
                   edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,          &
                   sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D,                  &
                   vdfg,maxmf,nupdraft,ktop_plume,                           &
                               tke_myj,el_myj,mixht,akhs,akms,elflx,ct,snow,sice,thz0,qz0,uz0,vz0,th10,q10,z0base,zntmyj,lowlyr,ivgtyp, &
                               thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten, &
                               icenter,jcenter,xcenter,ycenter,adaptmovetim,mvrec,nwritemv,ug,vg, &
                               gamk,kmw,ufw,vfw,u1b,v1b,cmz,csz,ce1z,ce2z,                  &
                               dum1,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p)
          ! end_write_restart
      use input
      use constants
      use lsnudge_module

      use netcdf
      use writeout_nc_module, only : restart_prelim,disp_err

      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes the writing of restart files
      !----------------------------------------------------------

      integer, intent(in) :: nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,trecs,trecw,arecs,arecw
      integer, intent(in) :: nwrite,nwritet,nwritea,nwriteh,nrst,nstatout
      integer, intent(in) :: num_soil_layers,nrad2d
      double precision, intent(in)    :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      real, intent(in) :: dt,dtlast
      integer, intent(in) :: ndt
      double precision, intent(in) :: adt,acfl,dbldt
      double precision, intent(in) :: mass1
      double precision, intent(in) :: mtime,stattim,taptim,rsttim,radtim,prcltim
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout), dimension(numq) :: asq,bsq
      character(len=3), intent(in), dimension(maxq) :: qname
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0,prs0,pi0,rho0,qv0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(in), dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, intent(in), dimension(ib:ie,jb:je) :: tsk,znt,ust,cd,ch,cq,u1,v1,s1,t1,xland,psfc,thflux,qvflux,prate,ustt,ut,vt,st
      real, intent(in), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(in), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,effc,effi,effs,effr,effg,effis
      real, intent(in), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(in), dimension(ni,nj,nrad2d) :: rad2d
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(in), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx,    &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta
      real, intent(in), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(in), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(in), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(in), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,3) :: ploc
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppx
      real, intent(in),    dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
      real, intent(in)   , dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: tsq,qsq,cov,sh3d,el_pbl,qc_bl,qi_bl,cldfra_bl, &
           qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,  &
           sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D
      real, intent(in)   , dimension(ibmynn:iemynn,jbmynn:jemynn) :: vdfg,maxmf
      integer, intent(in)   , dimension(ibmynn:iemynn,jbmynn:jemynn) :: nupdraft,ktop_plume
      real, intent(in),    dimension(ibmyj:iemyj,jbmyj:jemyj,kbmyj:kemyj) :: tke_myj,el_myj
      real, intent(in),    dimension(ibmyj:iemyj,jbmyj:jemyj) :: mixht,akhs,akms,elflx,ct,snow,sice,thz0,qz0,uz0,vz0,th10,q10,z0base,zntmyj
      integer, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: lowlyr,ivgtyp
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten
      real, intent(in   ) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(in   ) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      integer, intent(in   ) :: icenter,jcenter
      real, intent(in   ) :: xcenter,ycenter
      double precision, intent(in   ) :: adaptmovetim
      integer, intent(in   ) :: mvrec,nwritemv
      real, intent(in   ), dimension(kb:ke) :: ug,vg
      real, intent(in   ), dimension(kb:ke) :: gamk,kmw,ufw,vfw,u1b,v1b,cmz,csz,ce1z,ce2z
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
      integer, intent(inout), dimension(d3t) :: reqt
      integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

      character(len=maxstring) :: fname
      character(len=8) :: text1
      character(len=6) :: aname
      integer :: i,j,k,n,np,nvar,reqs,orecs,orecu,orecv,orecw,ndum
      integer :: ncid,time_index
      real, dimension(:), allocatable :: dumx,dumy
      integer :: nodemaster

      integer :: varid,ncstatus


!-----------------------------------------------------------------------

  IF( restart_format.eq.1 )THEN
    ! unformatted direct-access (grads) format:

  IF( restart_filetype.eq.1 )THEN

    !------------------
    ! one restart file (per stagger type):
    IF(myid.eq.nodemaster)THEN
      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_x.dat'

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=50,file=fname,form='unformatted',status='unknown')

      fname = 'cm1rst_s.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=51,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      orecs = 1

      fname = 'cm1rst_u.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=52,file=fname,form='unformatted',access='direct',recl=4*(nx+1)*ny)
      orecu = 1

      fname = 'cm1rst_v.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=53,file=fname,form='unformatted',access='direct',recl=4*nx*(ny+1))
      orecv = 1


      fname = 'cm1rst_w.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=54,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      orecw = 1

      if(dowr) write(outfile,*)
    ENDIF

  ELSEIF( restart_filetype.eq.2 )THEN

    !------------------
    ! one restart file (per restart time):
    IF(myid.eq.nodemaster)THEN
      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_x.dat'
      write(fname( 8:13),101) nrst
101   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=50,file=fname,form='unformatted',status='unknown')

      fname = 'cm1rst_XXXXXX_s.dat'
      write(fname( 8:13),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=51,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      orecs = 1

      fname = 'cm1rst_XXXXXX_u.dat'
      write(fname( 8:13),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=52,file=fname,form='unformatted',access='direct',recl=4*(nx+1)*ny)
      orecu = 1

      fname = 'cm1rst_XXXXXX_v.dat'
      write(fname( 8:13),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=53,file=fname,form='unformatted',access='direct',recl=4*nx*(ny+1))
      orecv = 1

      fname = 'cm1rst_XXXXXX_w.dat'
      write(fname( 8:13),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=54,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      orecw = 1

      if(dowr) write(outfile,*)
    ENDIF

  ELSEIF( restart_filetype.eq.3 )THEN

    !------------------
    ! one restart file per node (cm1r17 format):
    IF(myid.eq.nodemaster)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_YYYYYY.dat'

      write(fname( 8:13),102) mynode
      write(fname(15:20),102) nrst
102   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname
      if(dowr) write(outfile,*)

      open(unit=50,file=fname,form='unformatted',status='unknown')
    ENDIF
  ELSE
    stop 12388
  ENDIF


  ELSEIF( restart_format.eq.2 )THEN
    ! netcdf format:

    if( myid.eq.0 )then
      call     restart_prelim(nrst,ncid,mtime,xfref,yfref,zh,zf,sigma,sigmaf,  &
                              qname,num_soil_layers,nrad2d,dat2(1,1),dat2(1,2),dum1(ib,jb,kb),time_index)
    endif


  ELSE

    if( myid.eq.0 )then
      print *
      print *,'  unrecognized value for restart_format '
      print *
      print *,'      restart_format = ',restart_format
      print *
    endif

    call stopcm1

  ENDIF

!---------------------------------------------------------------
! metadata:

  IF( restart_format.eq.1 )THEN
    IF(myid.eq.0)THEN
      ! only processor 0 does this:
      write(50) nstep
      write(50) srec
      write(50) sirec
      write(50) urec
      write(50) vrec
      write(50) wrec
      write(50) nrec
      write(50) mrec
      write(50) prec
      write(50) trecs
      write(50) trecw
      write(50) arecs
      write(50) arecw
      write(50) mvrec
      write(50) nwrite
      write(50) nwritet
      write(50) nwritea
      write(50) nwritemv
      write(50) nrst
      write(50) nstatout
      write(50) ndt
      write(50) icenter
      write(50) jcenter
      write(50) output_format
      write(50) dt
      write(50) dtlast
      write(50) xcenter
      write(50) ycenter
      write(50) umove
      write(50) vmove
      write(50) adaptmovetim
      write(50) cflmax
      write(50) mtime
      write(50) stattim
      write(50) taptim
      write(50) rsttim
      write(50) radtim
      write(50) prcltim
      write(50) adt
      write(50) acfl
      write(50) dbldt
      write(50) mass1
      write(50) avgsfcu
      write(50) avgsfcv
      write(50) avgsfcs
      write(50) avgsfcsu
      write(50) avgsfcsv
      write(50) avgsfct
      write(50) avgsfcq
      write(50) avgsfcp
    ENDIF

  ELSEIF( restart_format.eq.2 )THEN
    IF(myid.eq.0)THEN

      call disp_err( nf90_inq_varid(ncid,"nstep",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nstep,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"srec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,srec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"sirec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,sirec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"urec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,urec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,vrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"wrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,wrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,prec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,trecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecw",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,trecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,arecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecw",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,arecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mvrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mvrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwrite",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwrite,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritet",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritet,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritea",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritea,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritemv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritemv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwriteh",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwriteh,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrst",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nrst,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nstatout",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nstatout,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ndt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,ndt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"icenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,icenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"jcenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,jcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"old_format",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,output_format,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dtlast",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"xcenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,xcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ycenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,ycenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umove",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,umove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmove",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,vmove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adaptmovetim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,adaptmovetim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"cflmax",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,cflmax,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mtime",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"stattim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,stattim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"taptim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,taptim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"rsttim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,rsttim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"radtim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,radtim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prcltim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,prcltim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,adt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"acfl",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,acfl,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dbldt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dbldt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mass1",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mass1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcu",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsu",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcsu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcsv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfct",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfct,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcq",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcq,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcp",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcp,(/time_index/)) , .true. )

    ENDIF

  ENDIF

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! budget variables:

    IF( myid.eq.0 )THEN

      IF( restart_format.eq.1 )THEN
        write(50) qbudget
        write(50) asq
        write(50) bsq

      ELSEIF( restart_format.eq.2 )THEN
        call disp_err( nf90_inq_varid(ncid,"qbudget",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,qbudget,(/1,time_index/),(/nbudget,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"asq",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,asq,(/1,time_index/),(/numq,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"bsq",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,bsq,(/1,time_index/),(/numq,1/)) , .true. )

      ENDIF

    ENDIF

!---------------------------------------------------------------

  IF( do_lsnudge .or. do_adapt_move )THEN

    if( myid.eq.0 )then
      write(50) ug
      write(50) vg
    endif

  ENDIF

!---------------------------------------------------------------
! standard 2D:

      n = 1
      call writer(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    if( nrain.eq.2 )then
      n = 2
      call writer(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain2   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    endif
      call writer(ni,nj,1,1,nx,ny,tsk(ib,jb),'tsk     ',           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)

!---------------------------------------------------------------
! standard 3D:

      call writer(ni,nj,1,nk,nx,ny,rho(ib,jb,1),'rho     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,prs(ib,jb,1),'prs     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),'ua      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2iu,d2ju,d3iu,d3ju)
      call writer(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),'va      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2iv,d2jv,d3iv,d3jv)
      call writer(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),'wa      ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),'ppi     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,tha(ib,jb,1),'tha     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,ppx(ib,jb,1),'ppx     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    if( psolver.eq.6 )then
      call writer(ni,nj,1,nk,nx,ny,phi1(ib,jb,1),'phi1    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,phi2(ib,jb,1),'phi2    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    endif
    IF(imoist.eq.1)THEN
    do n=1,numq
      text1 = '        '
      write(text1(1:3),156) qname(n)
156   format(a3)
      call writer(ni,nj,1,nk,nx,ny,qa(ib,jb,1,n),text1     ,       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    enddo
    ENDIF
    if(imoist.eq.1.and.eqtset.eq.2)then
      call writer(ni,nj,1,nk,nx,ny,qpten(ib,jb,1),'qpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qtten(ib,jb,1),'qtten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qvten(ib,jb,1),'qvten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qcten(ib,jb,1),'qcten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    endif
    if( cm1setup.eq.1 .and. iusetke )then
      call writer(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),'tkea    ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    endif

!---------------------------------------------------------------
!  radiation:

      if(radopt.ge.1)then
        call writer(ni,nj,1,nk,nx,ny,lwten(ib,jb,1),'lwten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,swten(ib,jb,1),'swten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radsw(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radsw   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rnflx(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'rnflx   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radswnet(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radswnet',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radlwin(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radlwin ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        do n=1,nrad2d
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rad2d(i,j,n)
        enddo
        enddo
        if( n.lt.10 )then
          text1 = 'radX    '
          write(text1(4:4),181) n
181       format(i1.1)
        elseif( n.lt.100 )then
          text1 = 'radXX   '
          write(text1(4:5),182) n
182       format(i2.2)
        elseif( n.lt.1000 )then
          text1 = 'radXXX  '
          write(text1(4:6),183) n
183       format(i3.3)
        else
          stop 11611
        endif
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        enddo
      endif

      if( radopt.ge.1 .and. ptype.eq.5 )then
        call writer(ni,nj,1,nk,nx,ny,effc(ib,jb,1),'effc    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effi(ib,jb,1),'effi    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effs(ib,jb,1),'effs    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effr(ib,jb,1),'effr    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effg(ib,jb,1),'effg    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effis(ib,jb,1),'effis   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      endif

!---------------------------------------------------------------
!  surface:
!     I don't know how many of these are really needed in restart
!     files, but let's include them all for now ... just to be safe

      if((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1))then
        !---- (1) ----!
      if(sfcmodel.ge.1)then
        call writer(ni,nj,1,1,nx,ny,ust(ib,jb),'ust     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,znt(ib,jb),'znt     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,cd(ib,jb),'cd      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ch(ib,jb),'ch      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,cq(ib,jb),'cq      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,u1(ib,jb),'u1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,v1(ib,jb),'v1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,s1(ib,jb),'s1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,t1(ib,jb),'t1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,u10(ib,jb),'u10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,v10(ib,jb),'v10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,s10(ib,jb),'s10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,xland(ib,jb),'xland   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,thflux(ib,jb),'thflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,qvflux(ib,jb),'qvflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psfc(ib,jb),'psfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      if( tbc.eq.3 )then
        call writer(ni,nj,1,1,nx,ny,ustt(ib,jb),'ustt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  ut(ib,jb),'ut      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  vt(ib,jb),'vt      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  st(ib,jb),'st      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      endif
      endif


      if(sfcmodel.ge.1)then
        !---- (2) ----!
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = lu_index(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lu_index',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = kpbl2d(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'kpbl2d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hfx(ib,jb),'hfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,qfx(ib,jb),'qfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hpbl(ib,jb),'hpbl    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,wspd(ib,jb),'wspd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psim(ib,jb),'psim    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psih(ib,jb),'psih    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,gz1oz0(ib,jb),'gz1oz0  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,br(ib,jb),'br      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CHS(ib,jb),'chs     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CHS2(ib,jb),'chs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CQS2(ib,jb),'cqs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CPMM(ib,jb),'cpmm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ZOL(ib,jb),'zol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,MAVAIL(ib,jb),'mavail  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,MOL(ib,jb),'mol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,RMOL(ib,jb),'rmol    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,REGIME(ib,jb),'regime  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,LH(ib,jb),'lh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,tmn(ib,jb),'tmn     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,FLHC(ib,jb),'flhc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,FLQC(ib,jb),'flqc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,QGH(ib,jb),'qgh     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CK(ib,jb),'ck      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CKA(ib,jb),'cka     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CDA(ib,jb),'cda     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,USTM(ib,jb),'ustm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,QSFC(ib,jb),'qsfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,T2(ib,jb),'t2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,Q2(ib,jb),'q2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,TH2(ib,jb),'th2     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,EMISS(ib,jb),'emiss   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,THC(ib,jb),'thc     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ALBD(ib,jb),'albd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,gsw(ib,jb),'gsw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,glw(ib,jb),'glw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,chklowq(ib,jb),'chklowq ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,capg(ib,jb),'capg    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,snowc(ib,jb),'snowc   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,fm(ib,jb),'fm      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,fh(ib,jb),'fh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,mznt(ib,jb),'mznt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,swspd(ib,jb),'swspd   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,wstar(ib,jb),'wstar   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,delta(ib,jb),'delta   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        do n=1,num_soil_layers
          if( n.lt.10 )then
            text1 = 'tslbX   '
            write(text1(5:5),171) n
171         format(i1.1)
          elseif( n.lt.100 )then
            text1 = 'tslbXX  '
            write(text1(5:6),172) n
172         format(i2.2)
          else
            stop 22122
          endif
          call writer(ni,nj,1,1,nx,ny,tslb(ib,jb,n),text1,             &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        enddo
      endif
      endif

      if(oceanmodel.eq.2)then
        call writer(ni,nj,1,1,nx,ny,tml(ib,jb),'tml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,t0ml(ib,jb),'t0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hml(ib,jb),'hml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,h0ml(ib,jb),'h0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,huml(ib,jb),'huml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hvml(ib,jb),'hvml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,tmoml(ib,jb),'tmoml   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      endif

!---------------------------------------------------------------

    IF( ipbl.eq.4 .or. ipbl.eq.5 .or. sfcmodel.eq.6 )THEN
      if(myid.eq.0) print *,'  writing mynn vars ... '
      call writer(ni,nj,1,nk,nx,ny,qke(ib,jb,1),'qke     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qke3d(ib,jb,1),'qke3d   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,tsq(ib,jb,1),'tsq     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qsq(ib,jb,1),'qsq     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,cov(ib,jb,1),'cov     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,sh3d(ib,jb,1),'sh3d    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,el_pbl(ib,jb,1),'el_pbl  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qc_bl(ib,jb,1),'qc_bl   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qi_bl(ib,jb,1),'qi_bl   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,cldfra_bl(ib,jb,1),'cldfra_b',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_a(ib,jb,1),'edmf_a  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_w(ib,jb,1),'edmf_w  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_qt(ib,jb,1),'edmf_qt ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_thl(ib,jb,1),'edmf_thl',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_ent(ib,jb,1),'edmf_ent',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_qc(ib,jb,1),'edmf_qc ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,sub_thl3d(ib,jb,1),'sub_thl3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,sub_sqv3d(ib,jb,1),'sub_sqv3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,det_thl3d(ib,jb,1),'det_thl3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,det_sqv3d(ib,jb,1),'det_sqv3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,thpten(ib,jb,1),'thpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qvpten(ib,jb,1),'qvpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qcpten(ib,jb,1),'qcpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qipten(ib,jb,1),'qipten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,upten(ib,jb,1),'upten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,vpten(ib,jb,1),'vpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qnipten(ib,jb,1),'qnipten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qncpten(ib,jb,1),'qncpten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,vdfg(ib,jb),'vdfg    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,maxmf(ib,jb),'maxmf   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = nupdraft(i,j)
        dum1(i,j,2) = ktop_plume(i,j)
      enddo
      enddo
      call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'nupdraft',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'ktop_plu',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)






      if(myid.eq.0) print *,'  ... done writing mynn vars '
    ENDIF

!---------------------------------------------------------------

    IF(ipbl.eq.6)THEN
      if(myid.eq.0) print *,'  writing myj vars ... '
      call writer(ni,nj,1,nk,nx,ny,tke_myj(ib,jb,1),'tke_myj ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,el_myj(ib,jb,1),'el_myj  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF(sfcmodel.eq.7)THEN
      call writer(ni,nj,1,1,nx,ny,mixht(ib,jb),'mixht   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,akhs(ib,jb),'akhs    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,akms(ib,jb),'akms    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,elflx(ib,jb),'elflx   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,ct(ib,jb),'ct      ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,snow(ib,jb),'snow    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sice(ib,jb),'sice    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,thz0(ib,jb),'thz0    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,qz0(ib,jb),'qz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,uz0(ib,jb),'uz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,vz0(ib,jb),'vz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,th10(ib,jb),'th10    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,q10(ib,jb),'q10     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,z0base(ib,jb),'z0base  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,zntmyj(ib,jb),'zntmyj  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = lowlyr(i,j)
        dum1(i,j,2) = ivgtyp(i,j)
      enddo
      enddo
      call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lowlyr  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'ivgtyp  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)






      if(myid.eq.0) print *,'  ... done writing myj vars '
    endif

!---------------------------------------------------------------
!  2-part turbulence param:

      if( myid.eq.0 )then
        write(50) gamk
        write(50) kmw
        write(50) ufw
        write(50) vfw
        write(50) u1b
        write(50) v1b
        write(50) cmz
        write(50) csz
        write(50) ce1z
        write(50) ce2z
      endif

!---------------------------------------------------------------
!  passive tracers:

      if(iptra.eq.1)then

        if(myid.eq.0)then
          if( restart_format.eq.1 )then
            write(50) npt

          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"npt",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,npt,(/time_index/)) , .true. )

          endif
        endif
        do n=1,npt
          if( n.lt.10 )then
            text1 = 'ptX     '
            write(text1(3:3),161) n
161         format(i1.1)
          elseif( n.lt.100 )then
            text1 = 'ptXX    '
            write(text1(3:4),162) n
162         format(i2.2)
          else
            stop 11512
          endif
          call writer(ni,nj,1,nk,nx,ny,pta(ib,jb,1,n),text1,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        enddo
      else
        if(myid.eq.0)then
          nvar = 0
          if( restart_format.eq.1 )then
            write(50) nvar

          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"npt",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,nvar,(/time_index/)) , .true. )

          endif
        endif
      endif

!---------------------------------------------------------------
!  parcels:

      if(iprcl.eq.1)then
        !-----------------
        ! with parcels:
        if(myid.eq.0)then
          if( restart_format.eq.1 )then
            write(50) nparcels

          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"numparcels",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,nparcels,(/time_index/)) , .true. )

          endif
        endif
        ! only write position info:
        if(myid.eq.0)then
          if( .not. terrain_flag )then
            DO np=1,nparcels
              ploc(np,1)=pdata(np,prx)
              ploc(np,2)=pdata(np,pry)
              ploc(np,3)=pdata(np,prz)
            ENDDO
          else
            DO np=1,nparcels
              ploc(np,1)=pdata(np,prx)
              ploc(np,2)=pdata(np,pry)
              ploc(np,3)=pdata(np,prsig)
            ENDDO
          endif
          if( restart_format.eq.1 )then
            write(50) ploc

          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
            n = 3
            call disp_err( nf90_put_var(ncid,varid,ploc,(/1,1,time_index/),(/nparcels,n,1/)) , .true. )

          endif
        endif
      else
        !-----------------
        ! without parcels:
        if(myid.eq.0)then
          nvar = 0
          if( restart_format.eq.1 )then
            write(50) nvar

          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"numparcels",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,nvar,(/time_index/)) , .true. )

          endif
        endif
        !-----------------
      endif

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  open bc:

      if(irbc.eq.4)then
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0)then
          ndum = ny
        else
          ndum = 1
        endif
        allocate( dumy(ndum) )
        !----------------------
      if( wbc.eq.2 )then
        aname = 'radbcw'

        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in writerbcwe, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif

        do k=1,nk
          call writerbcwe(radbcw,aname,ndum,dumy,ibw,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              write(50) dumy

            elseif( restart_format.eq.2 )then
              ncstatus = nf90_put_var(ncid,varid,dumy,(/1,k,time_index/),(/ny,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in writerbcwe, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif

            endif
          endif
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      if( ebc.eq.2 )then
        aname = 'radbce'

        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in writerbcwe, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif

        do k=1,nk
          call writerbcwe(radbce,aname,ndum,dumy,ibe,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              write(50) dumy

            elseif( restart_format.eq.2 )then
              ncstatus = nf90_put_var(ncid,varid,dumy,(/1,k,time_index/),(/ny,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in writerbcwe, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif

            endif
          endif
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        deallocate( dumy )
        if(myid.eq.0)then
          ndum = nx
        else
          ndum = 1
        endif
        allocate( dumx(ndum) )
        !----------------------
      if( sbc.eq.2 )then
        aname = 'radbcs'

        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in writerbcsn, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif

        do k=1,nk
          call writerbcsn(radbcs,aname,ndum,dumx,ibs,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              write(50) dumx

            elseif( restart_format.eq.2 )then
              ncstatus = nf90_put_var(ncid,varid,dumx,(/1,k,time_index/),(/nx,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in writerbcsn, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif

            endif
          endif
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      if( nbc.eq.2 )then
        aname = 'radbcn'

        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in writerbcsn, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif

        do k=1,nk
          call writerbcsn(radbcn,aname,ndum,dumx,ibn,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              write(50) dumx

            elseif( restart_format.eq.2 )then
              ncstatus = nf90_put_var(ncid,varid,dumx,(/1,k,time_index/),(/nx,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in writerbcsn, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif

            endif
          endif
        enddo
      endif
        !----------------------
        deallocate( dumx )
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      endif

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  150820:  optional variables


    IF( restart_file_theta )THEN
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k) = th0(i,j,k)+tha(i,j,k)
      enddo
      enddo
      enddo
      call writer(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'theta   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_u0 .or. do_lsnudge .or. do_adapt_move )THEN
      call writer(ni+1,nj,1,nk,nx+1,ny,u0(ib,jb,1),'u0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2iu,d2ju,d3iu,d3ju)
    ENDIF
    IF( restart_file_v0 .or. do_lsnudge .or. do_adapt_move )THEN
      call writer(ni,nj+1,1,nk,nx,ny+1,v0(ib,jb,1),'v0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2iv,d2jv,d3iv,d3jv)
    ENDIF
    IF( restart_file_dbz )THEN
      ! cm1r19:
      call writer(ni,nj,1,nk,nx,ny,qdiag(ibdq,jbdq,1,qd_dbz),'dbz     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    !-----
    IF( restart_file_th0 )THEN
      call writer(ni,nj,1,nk,nx,ny,th0(ib,jb,1),'th0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_prs0 )THEN
      call writer(ni,nj,1,nk,nx,ny,prs0(ib,jb,1),'prs0    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_pi0 )THEN
      call writer(ni,nj,1,nk,nx,ny,pi0(ib,jb,1),'pi0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_rho0 )THEN
      call writer(ni,nj,1,nk,nx,ny,rho0(ib,jb,1),'rho0    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_qv0 )THEN
      call writer(ni,nj,1,nk,nx,ny,qv0(ib,jb,1),'qv0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    !-----
    IF( restart_file_zs )THEN
      call writer(ni,nj,1,1,nx,ny,zs(ib,jb),'zs      ',            &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_zh )THEN
      call writer(ni,nj,1,nk,nx,ny,zh(ib,jb,1),'zhalf   ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_zf )THEN
      call writer(ni,nj,1,nk+1,nx,ny,zf(ib,jb,1),'zfull   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF

    IF( restart_file_diags )THEN
      if( td_diss.gt.0 )                                                 &
      call writer(ni,nj,1,nk,nx,ny,tdiag(ib,jb,1,td_diss),'dissheat',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      if( td_mp.gt.0 )                                                   &
      call writer(ni,nj,1,nk,nx,ny,tdiag(ib,jb,1,td_mp),'mptend  ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      if( qd_vtc.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtc),'vtc     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      if( qd_vtr.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtr),'vtr     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      if( qd_vts.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vts),'vts     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      if( qd_vtg.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtg),'vtg     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      if( qd_vti.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vti),'vti     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------

      ! cm1r19.7 fix !
      if( myid.eq.0 )then
      if( restart_format.eq.1 )then
        write(50) nwriteh

      elseif( restart_format.eq.2 )then
      call disp_err( nf90_inq_varid(ncid,"nwriteh",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwriteh,(/time_index/)) , .true. )

      endif
      endif




    IF( restart_format.eq.1 )THEN
      IF(myid.eq.0) close(unit=50)
      IF(myid.eq.nodemaster) close(unit=51)
      IF(myid.eq.nodemaster) close(unit=52)
      IF(myid.eq.nodemaster) close(unit=53)
      IF(myid.eq.nodemaster) close(unit=54)

    ELSEIF( restart_format.eq.2 )THEN
      if( myid.eq.0 )then
        call disp_err( nf90_close(ncid) , .true. )
      endif

    ENDIF



      return
      end subroutine write_restart


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine read_restart(nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,      &
                              trecs,trecw,arecs,arecw,                             &
                              nwrite,nwritet,nwritea,nwriteh,nrst,nstatout,        &
                              num_soil_layers,nrad2d,                              &
                              avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp, &
                              dt,dtlast,mtime,ndt,adt,acfl,dbldt,mass1,            &
                              stattim,taptim,rsttim,radtim,prcltim,                &
                              qbudget,asq,bsq,qname,                               &
                              xfref,yfref,zh,zf,sigma,sigmaf,zs,                   &
                              th0,prs0,pi0,rho0,qv0,u0,v0,                         &
                              rain,sws,svs,sps,srs,sgs,sus,shs,                    &
                              tsk,znt,ust,cd,ch,cq,u1,v1,s1,t1,thflux,qvflux,      &
                              prate,ustt,ut,vt,st,                                 &
                              radbcw,radbce,radbcs,radbcn,                         &
                              rho,prs,ua,va,wa,ppi,tha,qa,tkea,                    &
                              swten,lwten,radsw,rnflx,radswnet,radlwin,rad2d,      &
                              effc,effi,effs,effr,effg,effis,                      &
                              lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,xland,      &
                              hpbl,wspd,psim,psih,gz1oz0,br,                       &
                              CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                       &
                              MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                    &
                              CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,       &
                              gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta,tslb,    &
                              tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,               &
                              qpten,qtten,qvten,qcten,pta,pdata,ploc,ppx,          &
                              tdiag,qdiag,phi1,phi2,                               &
                   tsq,qsq,cov,sh3d,el_pbl,qc_bl,qi_bl,cldfra_bl,            &
                   qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,            &
                   edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,          &
                   sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D,                  &
                   vdfg,maxmf,nupdraft,ktop_plume,                           &
                              tke_myj,el_myj,mixht,akhs,akms,elflx,ct,snow,sice,thz0,qz0,uz0,vz0,th10,q10,z0base,zntmyj,lowlyr,ivgtyp, &
                              thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten, &
                              icenter,jcenter,xcenter,ycenter,adaptmovetim,mvrec,nwritemv,ug,vg, &
                              gamk,kmw,ufw,vfw,u1b,v1b,cmz,csz,ce1z,ce2z,                  &
                              dum1,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,restarted,restart_prcl)
          ! end_read_restart
      use input
      use constants
      use lsnudge_module
      use goddard_module, only : consat,consat2
      use lfoice_module, only : lfoice_init

      use netcdf
      use writeout_nc_module, only : disp_err

      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes the reading of restart files
      !----------------------------------------------------------

      integer, intent(inout) :: nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,trecs,trecw,arecs,arecw
      integer, intent(inout) :: nwrite,nwritet,nwritea,nwriteh,nrst,nstatout
      integer, intent(in) :: num_soil_layers,nrad2d
      double precision, intent(inout) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      real, intent(inout) :: dt,dtlast
      integer, intent(inout) :: ndt
      double precision, intent(inout) :: adt,acfl,dbldt
      double precision, intent(inout) :: mass1
      double precision, intent(inout) :: mtime,stattim,taptim,rsttim,radtim,prcltim
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout), dimension(numq) :: asq,bsq
      character(len=3), intent(in), dimension(maxq) :: qname
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0,prs0,pi0,rho0,qv0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,ust,cd,ch,cq,u1,v1,s1,t1,xland,psfc,thflux,qvflux,prate,ustt,ut,vt,st
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,effc,effi,effs,effr,effg,effis
      real, intent(inout), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(inout), dimension(ni,nj,nrad2d) :: rad2d
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(inout), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,3) :: ploc
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppx
      real, intent(inout), dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
      real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: tsq,qsq,cov,sh3d,el_pbl,qc_bl,qi_bl,cldfra_bl, &
           qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,  &
           sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D
      real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn) :: vdfg,maxmf
      integer, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn) :: nupdraft,ktop_plume
      real, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj,kbmyj:kemyj) :: tke_myj,el_myj
      real, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: mixht,akhs,akms,elflx,ct,snow,sice,thz0,qz0,uz0,vz0,th10,q10,z0base,zntmyj
      integer, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: lowlyr,ivgtyp
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten
      real, intent(in   ) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(in   ) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      integer, intent(inout) :: icenter,jcenter
      real, intent(inout) :: xcenter,ycenter
      double precision, intent(inout) :: adaptmovetim
      integer, intent(inout) :: mvrec,nwritemv
      real, intent(inout), dimension(kb:ke) :: ug,vg
      real, intent(inout), dimension(kb:ke) :: gamk,kmw,ufw,vfw,u1b,v1b,cmz,csz,ce1z,ce2z
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
      integer, intent(inout), dimension(d3t) :: reqt
      integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
      logical, intent(in) :: restarted
      logical, intent(inout) :: restart_prcl

      character(len=maxstring) fname
      character(len=8) :: text1
      character(len=6) :: aname
      integer :: i,j,k,n,np,nvar,nread,reqs,orecs,orecu,orecv,orecw,ndum
      integer :: ncid,time_index,old_format
      double precision, dimension(nbudget,0:numprocs-1) :: sbudget
      double precision, dimension(numq,0:numprocs-1) :: csq,dsq
      real, dimension(:,:), allocatable :: pfoo
      real, dimension(:), allocatable :: dumx,dumy
      integer :: nodemaster

      integer :: varid,ncstatus


!-----------------------------------------------------------------------

  IF( restart_format.eq.1 )THEN
    ! unformatted direct-access (grads) format:

  IF( restart_filetype.eq.1 )THEN

    !------------------
    ! one restart file (per stagger type):
    IF(myid.eq.nodemaster)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_x.dat'

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Reading from restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=50,file=fname,form='unformatted',status='old',err=778)

      fname = 'cm1rst_s.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=51,file=fname,form='unformatted',access='direct',recl=4*nx*ny,status='old',err=778)
      orecs = 1

      fname = 'cm1rst_u.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=52,file=fname,form='unformatted',access='direct',recl=4*(nx+1)*ny,status='old',err=778)
      orecu = 1

      fname = 'cm1rst_v.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=53,file=fname,form='unformatted',access='direct',recl=4*nx*(ny+1),status='old',err=778)
      orecv = 1

      fname = 'cm1rst_w.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=54,file=fname,form='unformatted',access='direct',recl=4*nx*ny,status='old',err=778)
      orecw = 1

      if(dowr) write(outfile,*)
    ENDIF

  ELSEIF( restart_filetype.eq.2 )THEN

    !------------------
    ! one restart file (per restart time):
    IF(myid.eq.nodemaster)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_x.dat'
      write(fname( 8:13),101) rstnum
101   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Reading from restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=50,file=fname,form='unformatted',status='old',err=778)

      fname = 'cm1rst_XXXXXX_s.dat'
      write(fname( 8:13),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=51,file=fname,form='unformatted',access='direct',recl=4*nx*ny,status='old',err=778)
      orecs = 1

      fname = 'cm1rst_XXXXXX_u.dat'
      write(fname( 8:13),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=52,file=fname,form='unformatted',access='direct',recl=4*(nx+1)*ny,status='old',err=778)
      orecu = 1

      fname = 'cm1rst_XXXXXX_v.dat'
      write(fname( 8:13),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=53,file=fname,form='unformatted',access='direct',recl=4*nx*(ny+1),status='old',err=778)
      orecv = 1

      fname = 'cm1rst_XXXXXX_w.dat'
      write(fname( 8:13),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=54,file=fname,form='unformatted',access='direct',recl=4*nx*ny,status='old',err=778)
      orecw = 1

      if(dowr) write(outfile,*)
    ENDIF

  ELSEIF( restart_filetype.eq.3 )THEN

    !------------------
    ! one restart file per node (cm1r17 format):
    IF(myid.eq.nodemaster)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_YYYYYY.dat'

      write(fname( 8:13),102) mynode
      write(fname(15:20),102) rstnum
102   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Reading from restart file!'
      if(dowr) write(outfile,*) '  fname=',fname
      if(dowr) write(outfile,*)

      open(unit=50,file=fname,form='unformatted',status='old')
    ENDIF
  ELSE
    stop 12389
  ENDIF


  ELSEIF( restart_format.eq.2 )THEN
    ! netcdf format:

    if( myid.eq.0 )then

      do i=1,maxstring
        string(i:i) = ' '
      enddo

    IF(     restart_filetype.eq.1 )THEN
      string = 'cm1rst.nc'
      time_index = rstnum
    ELSEIF( restart_filetype.eq.2 )THEN
      string = 'cm1rst_XXXXXX.nc'
      write(string(8:13),100) rstnum
100   format(i6.6)
      time_index = 1
    ENDIF
      if(myid.eq.0) print *,'  string = ',string

      call disp_err( nf90_open( path=string , mode=nf90_nowrite , ncid=ncid ) , .true. )

    endif


  ELSE

    if( myid.eq.0 )then
      print *
      print *,'  unrecognized value for restart_format '
      print *
      print *,'      restart_format = ',restart_format
      print *
    endif

    call stopcm1

  ENDIF

!---------------------------------------------------------------
! metadata:

  IF( restart_format.eq.1 )THEN
    IF(myid.eq.0)THEN
      ! only processor 0 has these variables:
      read(50) nstep
      read(50) srec
      read(50) sirec
      read(50) urec
      read(50) vrec
      read(50) wrec
      read(50) nrec
      read(50) mrec
      read(50) prec
      read(50) trecs
      read(50) trecw
      read(50) arecs
      read(50) arecw
      read(50) mvrec
      read(50) nwrite
      read(50) nwritet
      read(50) nwritea
      read(50) nwritemv
      read(50) nrst
      read(50) nstatout
      read(50) ndt
      read(50) icenter
      read(50) jcenter
      read(50) old_format
      read(50) dt
      read(50) dtlast
      read(50) xcenter
      read(50) ycenter
      read(50) umove
      read(50) vmove
      read(50) adaptmovetim
      read(50) cflmax
      read(50) mtime
      read(50) stattim
      read(50) taptim
      read(50) rsttim
      read(50) radtim
      read(50) prcltim
      read(50) adt
      read(50) acfl
      read(50) dbldt
      read(50) mass1
      read(50) avgsfcu
      read(50) avgsfcv
      read(50) avgsfcs
      read(50) avgsfcsu
      read(50) avgsfcsv
      read(50) avgsfct
      read(50) avgsfcq
      read(50) avgsfcp
    ENDIF

  ELSEIF( restart_format.eq.2 )THEN
    IF(myid.eq.0)THEN

      call disp_err( nf90_inq_varid(ncid,"nstep",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nstep,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"srec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,srec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"sirec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,sirec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"urec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,urec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,vrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"wrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,wrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,prec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,trecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecw",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,trecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,arecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecw",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,arecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mvrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mvrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwrite",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwrite,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritet",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritet,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritea",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritea,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritemv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritemv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwriteh",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwriteh,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrst",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nrst,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nstatout",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nstatout,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ndt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,ndt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"icenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,icenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"jcenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,jcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"old_format",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,old_format,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dtlast",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"xcenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,xcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ycenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,ycenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umove",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,umove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmove",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,vmove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adaptmovetim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,adaptmovetim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"cflmax",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,cflmax,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mtime",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"stattim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,stattim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"taptim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,taptim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"rsttim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,rsttim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"radtim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,radtim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prcltim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,prcltim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,adt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"acfl",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,acfl,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dbldt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dbldt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mass1",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mass1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcu",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsu",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcsu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcsv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfct",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfct,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcq",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcq,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcp",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcp,(/time_index/)) , .true. )

    ENDIF

  ENDIF



!---------------------------------------------------------------
! budget variables:

    IF( myid.eq.0 )THEN

      IF( restart_format.eq.1 )THEN
        read(50) qbudget
        read(50) asq
        read(50) bsq

      ELSEIF( restart_format.eq.2 )THEN
        call disp_err( nf90_inq_varid(ncid,"qbudget",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,qbudget,(/1,time_index/),(/nbudget,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"asq",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,asq,(/1,time_index/),(/numq,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"bsq",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,bsq,(/1,time_index/),(/numq,1/)) , .true. )

      ENDIF

    ELSE

      qbudget = 0.0
      asq = 0.0
      bsq = 0.0

    ENDIF

!---------------------------------------------------------------

  IF( do_lsnudge .or. do_adapt_move )THEN

    if( myid.eq.0 )then
      read(50) ug
      read(50) vg
    endif



  ENDIF

!---------------------------------------------------------------
! standard 2D:

      n = 1
      call  readr(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    if( nrain.eq.2 )then
      n = 2
      call  readr(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain2   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    endif
      call  readr(ni,nj,1,1,nx,ny,tsk(ib,jb),'tsk     ',           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)

!---------------------------------------------------------------
! standard 3D:

      call  readr(ni,nj,1,nk,nx,ny,rho(ib,jb,1),'rho     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,prs(ib,jb,1),'prs     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),'ua      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2iu,d2ju,d3iu,d3ju)
      call  readr(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),'va      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2iv,d2jv,d3iv,d3jv)
      call  readr(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),'wa      ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),'ppi     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,tha(ib,jb,1),'tha     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,ppx(ib,jb,1),'ppx     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    if( psolver.eq.6 )then
      call  readr(ni,nj,1,nk,nx,ny,phi1(ib,jb,1),'phi1    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,phi2(ib,jb,1),'phi2    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    endif
    IF(imoist.eq.1)THEN
    do n=1,numq
      text1 = '        '
      write(text1(1:3),156) qname(n)
156   format(a3)
      call  readr(ni,nj,1,nk,nx,ny,qa(ib,jb,1,n),text1     ,       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    enddo
    ENDIF
    if(imoist.eq.1.and.eqtset.eq.2)then
      call  readr(ni,nj,1,nk,nx,ny,qpten(ib,jb,1),'qpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qtten(ib,jb,1),'qtten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qvten(ib,jb,1),'qvten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qcten(ib,jb,1),'qcten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    endif
    if( cm1setup.eq.1 .and. iusetke )then
      call  readr(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),'tkea    ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    endif

!---------------------------------------------------------------
!  radiation:

      if(radopt.ge.1)then
        call  readr(ni,nj,1,nk,nx,ny,lwten(ib,jb,1),'lwten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,swten(ib,jb,1),'swten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radsw   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radsw(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'rnflx   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rnflx(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radswnet',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radswnet(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radlwin ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radlwin(i,j) = dum1(i,j,1)
        enddo
        enddo
        do n=1,nrad2d
        if( n.lt.10 )then
          text1 = 'radX    '
          write(text1(4:4),181) n
181       format(i1.1)
        elseif( n.lt.100 )then
          text1 = 'radXX   '
          write(text1(4:5),182) n
182       format(i2.2)
        elseif( n.lt.1000 )then
          text1 = 'radXXX  '
          write(text1(4:6),183) n
183       format(i3.3)
        else
          stop 11611
        endif
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rad2d(i,j,n) = dum1(i,j,1)
        enddo
        enddo
        enddo
      endif
      if( radopt.ge.1 .and. ptype.eq.5 )then
        call  readr(ni,nj,1,nk,nx,ny,effc(ib,jb,1),'effc    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effi(ib,jb,1),'effi    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effs(ib,jb,1),'effs    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effr(ib,jb,1),'effr    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effg(ib,jb,1),'effg    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effis(ib,jb,1),'effis   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      endif

!---------------------------------------------------------------
!  surface:
!     I don't know how many of these are really needed in restart
!     files, but let's include them all for now ... just to be safe

      if((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1))then
        !---- (1) ----!
      if(sfcmodel.ge.1)then
        call  readr(ni,nj,1,1,nx,ny,ust(ib,jb),'ust     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,znt(ib,jb),'znt     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,cd(ib,jb),'cd      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ch(ib,jb),'ch      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,cq(ib,jb),'cq      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,u1(ib,jb),'u1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,v1(ib,jb),'v1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,s1(ib,jb),'s1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,t1(ib,jb),'t1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,u10(ib,jb),'u10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,v10(ib,jb),'v10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,s10(ib,jb),'s10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,xland(ib,jb),'xland   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,thflux(ib,jb),'thflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,qvflux(ib,jb),'qvflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psfc(ib,jb),'psfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      if( tbc.eq.3 )then
        call  readr(ni,nj,1,1,nx,ny,ustt(ib,jb),'ustt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  ut(ib,jb),'ut      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  vt(ib,jb),'vt      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  st(ib,jb),'st      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      endif
      endif


      if(sfcmodel.ge.1)then
        !---- (2) ----!
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lu_index',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          lu_index(i,j) = nint(dum1(i,j,1))
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'kpbl2d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          kpbl2d(i,j) = nint(dum1(i,j,1))
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,hfx(ib,jb),'hfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,qfx(ib,jb),'qfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hpbl(ib,jb),'hpbl    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,wspd(ib,jb),'wspd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psim(ib,jb),'psim    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psih(ib,jb),'psih    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,gz1oz0(ib,jb),'gz1oz0  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,br(ib,jb),'br      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CHS(ib,jb),'chs     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CHS2(ib,jb),'chs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CQS2(ib,jb),'cqs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CPMM(ib,jb),'cpmm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ZOL(ib,jb),'zol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,MAVAIL(ib,jb),'mavail  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,MOL(ib,jb),'mol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,RMOL(ib,jb),'rmol    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,REGIME(ib,jb),'regime  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,LH(ib,jb),'lh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,tmn(ib,jb),'tmn     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,FLHC(ib,jb),'flhc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,FLQC(ib,jb),'flqc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,QGH(ib,jb),'qgh     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CK(ib,jb),'ck      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CKA(ib,jb),'cka     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CDA(ib,jb),'cda     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,USTM(ib,jb),'ustm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,QSFC(ib,jb),'qsfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,T2(ib,jb),'t2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,Q2(ib,jb),'q2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,TH2(ib,jb),'th2     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,EMISS(ib,jb),'emiss   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,THC(ib,jb),'thc     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ALBD(ib,jb),'albd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,gsw(ib,jb),'gsw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,glw(ib,jb),'glw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,chklowq(ib,jb),'chklowq ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,capg(ib,jb),'capg    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,snowc(ib,jb),'snowc   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,fm(ib,jb),'fm      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,fh(ib,jb),'fh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,mznt(ib,jb),'mznt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,swspd(ib,jb),'swspd   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,wstar(ib,jb),'wstar   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,delta(ib,jb),'delta   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        do n=1,num_soil_layers
          if( n.lt.10 )then
            text1 = 'tslbX   '
            write(text1(5:5),171) n
171         format(i1.1)
          elseif( n.lt.100 )then
            text1 = 'tslbXX  '
            write(text1(5:6),172) n
172         format(i2.2)
          else
            stop 22122
          endif
        call  readr(ni,nj,1,1,nx,ny,tslb(ib,jb,n),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        enddo
      endif
      endif

      if(oceanmodel.eq.2)then
        call  readr(ni,nj,1,1,nx,ny,tml(ib,jb),'tml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,t0ml(ib,jb),'t0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hml(ib,jb),'hml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,h0ml(ib,jb),'h0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,huml(ib,jb),'huml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hvml(ib,jb),'hvml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,tmoml(ib,jb),'tmoml   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      endif

!---------------------------------------------------------------

    IF( ipbl.eq.4 .or. ipbl.eq.5 .or. sfcmodel.eq.6 )THEN
      if(myid.eq.0) print *,'  reading mynn vars ... '
      call  readr(ni,nj,1,nk,nx,ny,qke(ib,jb,1),'qke     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qke3d(ib,jb,1),'qke3d   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,tsq(ib,jb,1),'tsq     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qsq(ib,jb,1),'qsq     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,cov(ib,jb,1),'cov     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,sh3d(ib,jb,1),'sh3d    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,el_pbl(ib,jb,1),'el_pbl  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qc_bl(ib,jb,1),'qc_bl   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qi_bl(ib,jb,1),'qi_bl   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,cldfra_bl(ib,jb,1),'cldfra_b',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_a(ib,jb,1),'edmf_a  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_w(ib,jb,1),'edmf_w  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_qt(ib,jb,1),'edmf_qt ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_thl(ib,jb,1),'edmf_thl',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_ent(ib,jb,1),'edmf_ent',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_qc(ib,jb,1),'edmf_qc ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,sub_thl3d(ib,jb,1),'sub_thl3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,sub_sqv3d(ib,jb,1),'sub_sqv3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,det_thl3d(ib,jb,1),'det_thl3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,det_sqv3d(ib,jb,1),'det_sqv3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,thpten(ib,jb,1),'thpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qvpten(ib,jb,1),'qvpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qcpten(ib,jb,1),'qcpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qipten(ib,jb,1),'qipten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,upten(ib,jb,1),'upten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,vpten(ib,jb,1),'vpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qnipten(ib,jb,1),'qnipten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qncpten(ib,jb,1),'qncpten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,vdfg(ib,jb),'vdfg    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,maxmf(ib,jb),'maxmf   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)






      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'nupdraft',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'ktop_plu',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        nupdraft(i,j) = nint(dum1(i,j,1))
        ktop_plume(i,j) = nint(dum1(i,j,2))
      enddo
      enddo
      if(myid.eq.0) print *,'  ... done reading mynn vars '
    ENDIF

!---------------------------------------------------------------

    IF(ipbl.eq.6)THEN
      if(myid.eq.0) print *,'  reading myj vars ... '
      call  readr(ni,nj,1,nk,nx,ny,tke_myj(ib,jb,1),'tke_myj ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,el_myj(ib,jb,1),'el_myj  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF(sfcmodel.eq.7)THEN
      call  readr(ni,nj,1,1,nx,ny,mixht(ib,jb),'mixht   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,akhs(ib,jb),'akhs    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,akms(ib,jb),'akms    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,elflx(ib,jb),'elflx   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,ct(ib,jb),'ct      ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,snow(ib,jb),'snow    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sice(ib,jb),'sice    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,thz0(ib,jb),'thz0    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,qz0(ib,jb),'qz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,uz0(ib,jb),'uz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,vz0(ib,jb),'vz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,th10(ib,jb),'th10    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,q10(ib,jb),'q10     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,z0base(ib,jb),'z0base  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,zntmyj(ib,jb),'zntmyj  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)






      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lowlyr  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'ivgtyp  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        lowlyr(i,j) = nint(dum1(i,j,1))
        ivgtyp(i,j) = nint(dum1(i,j,2))
      enddo
      enddo
      if(myid.eq.0) print *,'  ... done reading myj vars '
    endif

!---------------------------------------------------------------
!  2-part turbulence param:

      if( myid.eq.0 )then
         read(50) gamk
         read(50) kmw
         read(50) ufw
         read(50) vfw
         read(50) u1b
         read(50) v1b
         read(50) cmz
         read(50) csz
         read(50) ce1z
         read(50) ce2z
      endif


!---------------------------------------------------------------
!  passive tracers:

      if( myid.eq.0 )then
        if( restart_format.eq.1 )then
          read(50) nvar

        elseif( restart_format.eq.2 )then
          call disp_err( nf90_inq_varid(ncid,"npt",varid) , .true. )
          call disp_err( nf90_get_var(ncid,varid,nvar,(/time_index/)) , .true. )
          if( iptra.eq.0 ) nvar = 0

        endif
        print *,'  nvar_npt = ',nvar
      endif



      if( iptra.eq.1 .or. nvar.gt.0 )then
        if( nvar.gt.0 )then
          nread = 0
          if( iptra.eq.1 )then
            do n=1,min(nvar,npt)
              if( n.lt.10 )then
                text1 = 'ptX     '
                write(text1(3:3),161) n
161             format(i1.1)
              elseif( n.lt.100 )then
                text1 = 'ptXX    '
                write(text1(3:4),162) n
162             format(i2.2)
              else
                stop 11512
              endif
              call  readr(ni,nj,1,nk,nx,ny,pta(ib,jb,1,n),text1,           &
                          ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                          ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                          dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
              nread = nread+1
            enddo
          endif
          if( nread .lt. nvar )then
            ! need to read more data ....
            do n=nread+1,nvar
              if( n.lt.10 )then
                text1 = 'ptX     '
                write(text1(3:3),161) n
              elseif( n.lt.100 )then
                text1 = 'ptXX    '
                write(text1(3:4),162) n
              else
                stop 11513
              endif
              call  readr(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),text1,            &
                          ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                          ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                          dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
            enddo
          endif
        else
          if( myid.eq.0 ) print *
          if( myid.eq.0 ) print *,'  Note:  no passive tracer data in the restart file '
          if( myid.eq.0 ) print *
        endif
      endif

!---------------------------------------------------------------
!  parcels:

      if( myid.eq.0 )then
        if( restart_format.eq.1 )then
          read(50) nvar

        elseif( restart_format.eq.2 )then
          call disp_err( nf90_inq_varid(ncid,"numparcels",varid) , .true. )
          call disp_err( nf90_get_var(ncid,varid,nvar,(/time_index/)) , .true. )
          if( iprcl.eq.0 ) nvar = 0

        endif
        print *,'  nvar_parcels = ',nvar
      endif



      if( iprcl.eq.1 .or. nvar.gt.0 )then
        if( nvar.gt.0 )then
          ! only read position info:
          if( myid.eq.0 )then
            IF( nvar.eq.nparcels )THEN
              ! easy:  restart file matches current config
              if( restart_format.eq.1 )then
                read(50) ploc

              elseif( restart_format.eq.2 )then
                call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
                n = 3
                call disp_err( nf90_get_var(ncid,varid,ploc,(/1,1,time_index/),(/nparcels,n,1/)) , .true. )

              endif
            ELSE
              ! annoying:  restart file has different nparcels than current config
              IF( iprcl.eq.1 )THEN
                if( .not. terrain_flag )then
                  do np=1,nparcels
                    ploc(np,1) = pdata(np,prx)
                    ploc(np,2) = pdata(np,pry)
                    ploc(np,3) = pdata(np,prz)
                  enddo
                else
                  do np=1,nparcels
                    ploc(np,1) = pdata(np,prx)
                    ploc(np,2) = pdata(np,pry)
                    ploc(np,3) = pdata(np,prsig)
                  enddo
                endif
              ENDIF
              if( myid.eq.0 ) print *,'  start pfoo ' 
              allocate( pfoo(nvar,3) )
              if( restart_format.eq.1 )then
                read(50) pfoo

              elseif( restart_format.eq.2 )then
                call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
                n = 3
                call disp_err( nf90_get_var(ncid,varid,pfoo,(/1,1,time_index/),(/nvar,n,1/)) , .true. )

              endif
              IF( iprcl.eq.1 )THEN
                do n=1,3
                do np=1,min(nvar,nparcels)
                  ploc(np,n) = pfoo(np,n)
                enddo
                enddo
              ENDIF
              deallocate( pfoo )
              if( myid.eq.0 ) print *,'  end pfoo ' 
            ENDIF
          endif
          IF( iprcl.eq.1 )THEN

            if( .not. terrain_flag )then
              DO np=1,nparcels
                pdata(np,prx)=ploc(np,1)
                pdata(np,pry)=ploc(np,2)
                pdata(np,prz)=ploc(np,3)
              ENDDO
            else
              DO np=1,nparcels
                pdata(np,prx)=ploc(np,1)
                pdata(np,pry)=ploc(np,2)
                pdata(np,prsig)=ploc(np,3)
              ENDDO
            endif
            restart_prcl = .true.
          ENDIF
        else
          if( myid.eq.0 ) print *
          if( myid.eq.0 ) print *,'  Note:  no parcel data in the restart file '
          if( myid.eq.0 ) print *
          restart_prcl = .false.
        endif
      endif

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  open bc:

      if(irbc.eq.4)then
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0)then
          ndum = ny
        else
          ndum = 1
        endif
        allocate( dumy(ndum) )
        !----------------------
      if( wbc.eq.2 )then
        aname = 'radbcw'

        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in readrbcwe, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif

        do k=1,nk
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              read(50) dumy

            elseif( restart_format.eq.2 )then
              ncstatus = nf90_get_var(ncid,varid,dumy,(/1,k,time_index/),(/ny,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in readrbcwe, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif

            endif
          endif
          call readrbcwe(radbcw,aname,ndum,dumy,ibw,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      if( ebc.eq.2 )then
        aname = 'radbce'

        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in readrbcwe, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif

        do k=1,nk
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              read(50) dumy

            elseif( restart_format.eq.2 )then
              ncstatus = nf90_get_var(ncid,varid,dumy,(/1,k,time_index/),(/ny,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in readrbcwe, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif

            endif
          endif
          call readrbcwe(radbce,aname,ndum,dumy,ibe,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        deallocate( dumy )
        if(myid.eq.0)then
          ndum = nx
        else
          ndum = 1
        endif
        allocate( dumx(ndum) )
        !----------------------
      if( sbc.eq.2 )then
        aname = 'radbcs'

        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in readrbcsn, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif

        do k=1,nk
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              read(50) dumx

            elseif( restart_format.eq.2 )then
              ncstatus = nf90_get_var(ncid,varid,dumx,(/1,k,time_index/),(/nx,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in readrbcsn, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif

            endif
          endif
          call readrbcsn(radbcs,aname,ndum,dumx,ibs,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      if( nbc.eq.2 )then
        aname = 'radbcn'

        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in readrbcsn, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif

        do k=1,nk
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              read(50) dumx

            elseif( restart_format.eq.2 )then
              ncstatus = nf90_get_var(ncid,varid,dumx,(/1,k,time_index/),(/nx,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in readrbcsn, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif

            endif
          endif
          call readrbcsn(radbcn,aname,ndum,dumx,ibn,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
        enddo
      endif
        !----------------------
        deallocate( dumx )
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      endif

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  151001:  use theta  (over-rides perturbation value that was read-in above)


    IF( restart_use_theta )THEN
      call  readr(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'theta   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2is,d2js,d3is,d3js)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        tha(i,j,k) = dum1(i,j,k)-th0(i,j,k)
      enddo
      enddo
      enddo
    ENDIF

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------

    IF( do_lsnudge .or. do_adapt_move )THEN
      call  readr(ni+1,nj,1,nk,nx+1,ny,u0(ib,jb,1),'u0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2iu,d2ju,d3iu,d3ju)
      call  readr(ni,nj+1,1,nk,nx,ny+1,v0(ib,jb,1),'v0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2iv,d2jv,d3iv,d3jv)
    ENDIF

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------

      ! cm1r19.7 fix !
      if( myid.eq.0 )then
      if( restart_format.eq.1 )then
        read(50) nwriteh

      elseif( restart_format.eq.2 )then
      call disp_err( nf90_inq_varid(ncid,"nwriteh",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwriteh,(/time_index/)) , .true. )

      endif
      endif


    IF( restart_format.eq.1 )THEN
      IF(myid.eq.0) close(unit=50)
      IF(myid.eq.nodemaster) close(unit=51)
      IF(myid.eq.nodemaster) close(unit=52)
      IF(myid.eq.nodemaster) close(unit=53)
      IF(myid.eq.nodemaster) close(unit=54)

    ELSEIF( restart_format.eq.2 )THEN
      if( myid.eq.0 )then
        call disp_err( nf90_close(ncid) , .true. )
      endif

    ENDIF

    if( restarted ) nrst = nrst+1

!---------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  From restart file: '
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   mtime   = ',mtime
      if(dowr) write(outfile,*) '   stattim = ',stattim
      if(dowr) write(outfile,*) '   taptim  = ',taptim
      if(dowr) write(outfile,*) '   rsttim  = ',rsttim
      if(dowr) write(outfile,*) '   radtim  = ',radtim
      if(dowr) write(outfile,*) '   prcltim = ',prcltim
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   nstep   = ',nstep
      if(dowr) write(outfile,*) '   srec    = ',srec
      if(dowr) write(outfile,*) '   sirec   = ',sirec
      if(dowr) write(outfile,*) '   urec    = ',urec
      if(dowr) write(outfile,*) '   vrec    = ',vrec
      if(dowr) write(outfile,*) '   wrec    = ',wrec
      if(dowr) write(outfile,*) '   nrec    = ',nrec
      if(dowr) write(outfile,*) '   mrec    = ',mrec
      if(dowr) write(outfile,*) '   prec    = ',prec
      if(dowr) write(outfile,*) '   nwrite  = ',nwrite
      if(dowr) write(outfile,*) '   nrst    = ',nrst
      if(dowr) write(outfile,*) '   nstatout= ',nstatout
      if( dodomaindiag )then
      if(dowr) write(outfile,*) '   trecs   = ',trecs
      if(dowr) write(outfile,*) '   trecw   = ',trecw
      if(dowr) write(outfile,*) '   nwritet = ',nwritet
      endif
      if( doazimavg )then
      if(dowr) write(outfile,*) '   arecs   = ',arecs
      if(dowr) write(outfile,*) '   arecw   = ',arecw
      if(dowr) write(outfile,*) '   nwritea = ',nwritea
      endif
      if( dohifrq )then
      if(dowr) write(outfile,*) '   nwriteh = ',nwriteh
      endif
      if( do_adapt_move )then
      if(dowr) write(outfile,*) '   nwritemv = ',nwritemv
      if(dowr) write(outfile,*) '   mvrec    = ',mvrec
      endif
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   mass1 = ',mass1
      if(dowr) write(outfile,*)

!---------

      if( adapt_dt.eq.0 )then
        dt = dtl
        dbldt = dtl
      endif

      ! this is needed for stats files:
      nrec=nrec-1
      nstatout=nstatout-1

      IF( output_format .ne. old_format )THEN
        srec = 1
        sirec = 1
        urec = 1
        vrec = 1
        wrec = 1
        nrec = 1
        mrec = 1
        nwrite = 1
        prec = 1
      ENDIF

!---------

        IF( (imoist.eq.1).and.(ptype.eq.2) )then
          if(timestats.ge.1) time_misc=time_misc+mytime()
          call consat2(dt)
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF
        IF( (imoist.eq.1).and.(ptype.eq.4) )then
          if(timestats.ge.1) time_misc=time_misc+mytime()
          call lfoice_init(dt)
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF

!---------



      return

778   print *,'  error opening restart file '
      print *,'    ... stopping cm1 ... '
      call stopcm1

      end subroutine read_restart


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writerbcwe(radbc,aname,ndum,dumy,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)

      implicit none

      integer, intent(in) :: ndum,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,k
      real, intent(in), dimension(jb:je,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum) :: dumy
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(numprocs) :: myj1p

      integer :: j,j1,j2



      do j=1,ny
        dumy(j) = radbc(j,k)
      enddo


      end subroutine writerbcwe


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writerbcsn(radbc,aname,ndum,dumx,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)

      implicit none

      integer, intent(in) :: ndum,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,k
      real, intent(in), dimension(ib:ie,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum) :: dumx
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(numprocs) :: myi1p

      integer :: i,i1,i2



      do i=1,nx
        dumx(i) = radbc(i,k)
      enddo


      end subroutine writerbcsn


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine  readrbcwe(radbc,aname,ndum,dumy,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)

      implicit none

      integer, intent(in) :: ndum,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,k
      real, intent(inout), dimension(jb:je,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum) :: dumy
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(numprocs) :: myj1p

      integer :: j,j1,j2



      do j=1,ny
        radbc(j,k) = dumy(j)
      enddo


      end subroutine  readrbcwe


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine  readrbcsn(radbc,aname,ndum,dumx,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)

      implicit none

      integer, intent(in) :: ndum,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,k
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum) :: dumx
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(numprocs) :: myi1p

      integer :: i,i1,i2



      do i=1,nx
        radbc(i,k) = dumx(i)
      enddo


      end subroutine  readrbcsn


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writer(numi,numj,numk1,numk2,nxr,nyr,var,aname,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1,dat2,dat3,reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2i,d2j,d3i,d3j)

    use netcdf

    implicit none

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) and does the actual writing of restart files.
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(in   ), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=8), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: time_index,restart_format,restart_filetype
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in) :: mynode,nodemaster,nodes,nfile
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

    integer :: i,j,k,msk,nitmp,njtmp

    integer :: varid,status


!-------------------------------------------------------------------------------

    rf1:  IF( restart_filetype.eq.1 .or. restart_filetype.eq.2 )THEN

    if(myid.eq.0) print *,aname

    msk = 0



    kloop:  DO k=numk1,numk2


      !-------------------- non-MPI section --------------------!
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,numj
      do i=1,numi
        dat2(i,j)=var(i,j,k)
      enddo
      enddo

      if( restart_format.eq.2 )then
        status = nf90_inq_varid(ncid,aname,varid)
        if(status.ne.nf90_noerr)then
          print *,'  Error1 in writer, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
      endif


        ! WRITE DATA:
        IF( myid.eq.msk )THEN
          !---   write data   ------------------!
          IF( restart_format.eq.1 )THEN
            write(nfile,rec=orec) ((dat2(i,j),i=1,nxr),j=1,nyr)

          ELSEIF( restart_format.eq.2 )THEN
            ! ----- netcdf format -----
            if(numk1.eq.numk2)then
              status = nf90_put_var(ncid,varid,dat2,(/1,1,time_index/),(/d2i,d2j,1/))
            else
              status = nf90_put_var(ncid,varid,dat2,(/1,1,k,time_index/),(/d2i,d2j,1,1/))
            endif
            if(status.ne.nf90_noerr)then
              print *,'  Error2 in writer, aname = ',aname
              print *,nf90_strerror(status)
              call stopcm1
            endif

          ENDIF
        ENDIF

      !---  prepare for next level   -------!
      IF( restart_format.eq.1 )THEN
        orec = orec+1
!!!#ifdef MPI
!!!        msk = msk+ppnode
!!!        if( msk.ge.numprocs ) msk = msk-numprocs
!!!#endif
      ENDIF

      !---  done with this level   ---------!
    ENDDO  kloop

    ENDIF  rf1

!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
!ccccc  done  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-------------------------------------------------------------------------------

!!!#ifdef MPI
!!!    ! helps with memory:
!!!    call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!    !----------------- end MPI section -----------------!
!!!#endif

    return
    end subroutine writer


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine  readr(numi,numj,numk1,numk2,nxr,nyr,var,aname,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1,dat2,dat3,reqt,ppnode,d3n,d3t,mynode,nodemaster,nodes,d2i,d2j,d3i,d3j)

    use netcdf

    implicit none

    !-------------------------------------------------------------------
    ! This subroutine reads restart files and then passes data 
    ! to other processors if this is a MPI run. 
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=8), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: time_index,restart_format,restart_filetype
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in) :: mynode,nodemaster,nodes,nfile
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

    integer :: i,j,k,msk,nitmp,njtmp

    integer :: varid,status


!-------------------------------------------------------------------------------

    rf1:  IF( restart_filetype.eq.1 .or. restart_filetype.eq.2 )THEN

    if(myid.eq.0) print *,aname

    msk = 0


    if( myid.eq.0 )then
      if( restart_format.eq.2 )then
        status = nf90_inq_varid(ncid,aname,varid)
        if(status.ne.nf90_noerr)then
          print *,'  Error1 in  readr, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
      endif
    endif


    kloop:  DO k=numk1,numk2



          IF( restart_format.eq.1 )THEN
            read(nfile,rec=orec) ((dat2(i,j),i=1,nxr),j=1,nyr)

          ELSEIF( restart_format.eq.2 )THEN
            ! ----- netcdf format -----
            if(numk1.eq.numk2)then
              status = nf90_get_var(ncid,varid,dat2,(/1,1,time_index/),(/d2i,d2j,1/))
            else
              status = nf90_get_var(ncid,varid,dat2,(/1,1,k,time_index/),(/d2i,d2j,1,1/))
            endif
            if(status.ne.nf90_noerr)then
              print *,'  Error2 in  readr, aname = ',aname
              print *,nf90_strerror(status)
              call stopcm1
            endif

          ENDIF


!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,numj
          do i=1,numi
            var(i,j,k)=dat2(i,j)
          enddo
          enddo

      !---  prepare for next level   -------!
      IF( restart_format.eq.1 )THEN
        orec = orec+1
!!!#ifdef MPI
!!!        msk = msk+ppnode
!!!        if( msk.ge.numprocs ) msk = msk-numprocs
!!!#endif
      ENDIF

      !---  done with this level   ---------!
    ENDDO  kloop

    ENDIF  rf1

!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
!ccccc  done  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-------------------------------------------------------------------------------

!!!#ifdef MPI
!!!    ! helps with memory:
!!!    call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!    !----------------- end MPI section -----------------!
!!!#endif

    return
    end subroutine  readr




  END MODULE restart_module

