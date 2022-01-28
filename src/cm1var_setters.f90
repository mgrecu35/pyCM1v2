subroutine set_u3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,u3d_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: u3d_in(ibmo:iemo+1,jbmo:jemo,kbmo:kemo)
u3d=u3d_in
end subroutine set_u3d
!------------------
subroutine set_v3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,v3d_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: v3d_in(ibmo:iemo,jbmo:jemo+1,kbmo:kemo)
v3d=v3d_in
end subroutine set_v3d
!------------------
subroutine set_w3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,w3d_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: w3d_in(ibmo:iemo,jbmo:jemo,kbmo:kemo+1)
w3d=w3d_in
end subroutine set_w3d
!------------------
subroutine set_ua(ibmo,iemo,jbmo,jemo,kbmo,kemo,ua_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: ua_in(ibmo:iemo+1,jbmo:jemo,kbmo:kemo)
ua=ua_in
end subroutine set_ua
!------------------
subroutine set_va(ibmo,iemo,jbmo,jemo,kbmo,kemo,va_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: va_in(ibmo:iemo,jbmo:jemo+1,kbmo:kemo)
va=va_in
end subroutine set_va
!------------------
subroutine set_wa(ibmo,iemo,jbmo,jemo,kbmo,kemo,wa_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: wa_in(ibmo:iemo,jbmo:jemo,kbmo:kemo+1)
wa=wa_in
end subroutine set_wa
!------------------
subroutine set_ppi(ibmo,iemo,jbmo,jemo,kbmo,kemo,ppi_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: ppi_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
ppi=ppi_in
end subroutine set_ppi
!------------------
subroutine set_pp3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,pp3d_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: pp3d_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
pp3d=pp3d_in
end subroutine set_pp3d
!------------------
!---------Getters--------
subroutine get_u3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,u3d_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: u3d_out(ibmo:iemo+1,jbmo:jemo,kbmo:kemo)
u3d_out=u3d
end subroutine get_u3d
!------------------
subroutine get_v3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,v3d_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: v3d_out(ibmo:iemo,jbmo:jemo+1,kbmo:kemo)
v3d_out=v3d
end subroutine get_v3d
!------------------
subroutine get_w3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,w3d_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: w3d_out(ibmo:iemo,jbmo:jemo,kbmo:kemo+1)
w3d_out=w3d
end subroutine get_w3d
!------------------
subroutine get_ua(ibmo,iemo,jbmo,jemo,kbmo,kemo,ua_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: ua_out(ibmo:iemo+1,jbmo:jemo,kbmo:kemo)
ua_out=ua
end subroutine get_ua
!------------------
subroutine get_va(ibmo,iemo,jbmo,jemo,kbmo,kemo,va_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: va_out(ibmo:iemo,jbmo:jemo+1,kbmo:kemo)
va_out=va
end subroutine get_va
!------------------
subroutine get_wa(ibmo,iemo,jbmo,jemo,kbmo,kemo,wa_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: wa_out(ibmo:iemo,jbmo:jemo,kbmo:kemo+1)
wa_out=wa
end subroutine get_wa
!------------------
subroutine get_ppi(ibmo,iemo,jbmo,jemo,kbmo,kemo,ppi_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: ppi_out(ibmo:iemo,jbmo:jemo,kbmo:kemo)
ppi_out=ppi
end subroutine get_ppi
!------------------
subroutine get_pp3d(ibmo,iemo,jbmo,jemo,kbmo,kemo,pp3d_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: pp3d_out(ibmo:iemo,jbmo:jemo,kbmo:kemo)
pp3d_out=pp3d
end subroutine get_pp3d
!------------------
