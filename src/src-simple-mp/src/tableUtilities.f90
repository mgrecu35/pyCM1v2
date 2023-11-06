subroutine getqvs(t,press,qvs)
    use tabulated_mp
    implicit none
    integer :: it,ip
    real :: ft,fp,tmin,tmax,pmin,pmax,dtemp,dpres
    real :: t,press
    real :: qvs
    tmin=tdef(1)
    tmax=tdef(nt)
    pmin=pdef(1)
    pmax=pdef(np)
    if(t.ge.tmin.and.t.le.tmax.and.press.ge.pmin.and.press.le.pmax) then
        dtemp=tdef(2)-tdef(1)
        dpres=pdef(2)-pdef(1)
        it=int((t-tmin)/dtemp)+1
        ip=int((press-pmin)/dpres)+1
        ft=(t-tdef(it))/(tdef(it+1)-tdef(it))
        fp=(press-pdef(ip))/(pdef(ip+1)-pdef(ip))
        
        qvs=(1.-ft)*(1.-fp)*qvsat_map(it,ip)+ft*(1.-fp)*qvsat_map(it+1,ip)+ &
            (1.-ft)*fp*qvsat_map(it,ip+1)+ft*fp*qvsat_map(it+1,ip+1)
        !print*,t,press,it,ip,ft,fp,qvsat_map(it,ip),qvs
    else
        qvs=999.
    end if
end subroutine getqvs

subroutine getqvi(t,press,qvi)
    use tabulated_mp
    implicit none
    integer :: it,ip
    real :: ft,fp,tmin,tmax,pmin,pmax,dtemp,dpres
    real :: t,press
    real :: qvi
    tmin=tdef(1)
    tmax=tdef(nt)
    pmin=pdef(1)
    pmax=pdef(np)
    if(t.ge.tmin.and.t.le.tmax.and.press.ge.pmin.and.press.le.pmax) then
        dtemp=tdef(2)-tdef(1)
        dpres=pdef(2)-pdef(1)
        it=int((t-tmin)/dtemp)+1
        ip=int((press-pmin)/dpres)+1
        ft=(t-tdef(it))/(tdef(it+1)-tdef(it))
        fp=(press-pdef(ip))/(pdef(ip+1)-pdef(ip))
        
        qvi=(1.-ft)*(1.-fp)*qvsat_ice_map(it,ip)+ft*(1.-fp)*qvsat_ice_map(it+1,ip)+ &
            (1.-ft)*fp*qvsat_ice_map(it,ip+1)+ft*fp*qvsat_ice_map(it+1,ip+1)
        !print*,t,press,it,ip,ft,fp,qvsat_map(it,ip),qvs
    else
        qvi=999.
    end if
end subroutine getqvi

subroutine get_condensation(t,press,supersat,fract_tables)
    use tabulated_mp
    implicit none
    integer :: it,ip, ifract
    real :: ft,fp,tmin,tmax,pmin,pmax,dtemp,dpres
    real :: t,press
    real :: supersat, fract_tables, fract1, fract2, ffract
    tmin=tdef(1)
    tmax=tdef(nt)
    pmin=pdef(1)
    pmax=pdef(np)
    ifract=(supersat/0.01)
    if (ifract>20) ifract=20
    if (ifract<=0) then
        fract_tables=0.
        return
    end if
    if(t.ge.tmin.and.t.le.tmax.and.press.ge.pmin.and.press.le.pmax) then
        dtemp=tdef(2)-tdef(1)
        dpres=pdef(2)-pdef(1)
        it=int((t-tmin)/dtemp)+1
        ip=int((press-pmin)/dpres)+1
        ft=(t-tdef(it))/(tdef(it+1)-tdef(it))
        fp=(press-pdef(ip))/(pdef(ip+1)-pdef(ip))
        
        fract1=(1.-ft)*(1.-fp)*fract_map(it,ip,ifract)+ft*(1.-fp)*fract_map(it+1,ip,ifract)+ &
            (1.-ft)*fp*fract_map(it,ip+1,ifract)+ft*fp*fract_map(it+1,ip+1,ifract)
        if(ifract<20) then
            fract2=(1.-ft)*(1.-fp)*fract_map(it,ip,ifract+1)+ft*(1.-fp)*fract_map(it+1,ip,ifract+1)+ &
            (1.-ft)*fp*fract_map(it,ip+1,ifract+1)+ft*fp*fract_map(it+1,ip+1,ifract+1)
            ffract=(supersat-(ifract-1)*0.01)/0.01
            fract_tables=(fract1*(1.-ffract)+fract2*ffract)*0.01
        else
            fract_tables=fract1*0.01
        end if
        !print*,t,press,it,ip,ft,fp,qvsat_map(it,ip),qvs
    else
        fract_tables=0.
    end if
end subroutine get_condensation

subroutine get_sublimation(t,press,supersat,fract_tables)
    use tabulated_mp
    implicit none
    integer :: it,ip, ifract
    real :: ft,fp,tmin,tmax,pmin,pmax,dtemp,dpres
    real :: t,press
    real :: supersat, fract_tables, fract1, fract2, ffract
    tmin=tdef(1)
    tmax=tdef(nt)
    pmin=pdef(1)
    pmax=pdef(np)
    ifract=(supersat/0.01)
    if (ifract>20) ifract=20
    if (ifract<=0) then
        fract_tables=0.
        return
    end if
    if(t.ge.tmin.and.t.le.tmax.and.press.ge.pmin.and.press.le.pmax) then
        dtemp=tdef(2)-tdef(1)
        dpres=pdef(2)-pdef(1)
        it=int((t-tmin)/dtemp)+1
        ip=int((press-pmin)/dpres)+1
        ft=(t-tdef(it))/(tdef(it+1)-tdef(it))
        fp=(press-pdef(ip))/(pdef(ip+1)-pdef(ip))
        
        fract1=(1.-ft)*(1.-fp)*fract_ice_map(it,ip,ifract)+ft*(1.-fp)*fract_ice_map(it+1,ip,ifract)+ &
            (1.-ft)*fp*fract_ice_map(it,ip+1,ifract)+ft*fp*fract_ice_map(it+1,ip+1,ifract)
        if(ifract<20) then
            fract2=(1.-ft)*(1.-fp)*fract_ice_map(it,ip,ifract+1)+ft*(1.-fp)*fract_ice_map(it+1,ip,ifract+1)+ &
            (1.-ft)*fp*fract_ice_map(it,ip+1,ifract+1)+ft*fp*fract_ice_map(it+1,ip+1,ifract+1)
            ffract=(supersat-(ifract-1)*0.01)/0.01
            fract_tables=(fract1*(1.-ffract)+fract2*ffract)*0.01
        else
            fract_tables=fract1*0.01
        end if
        !print*,t,press,it,ip,ft,fp,qvsat_map(it,ip),qvs
    else
        fract_tables=0.
    end if
end subroutine get_sublimation