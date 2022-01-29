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
subroutine set_prs(ibmo,iemo,jbmo,jemo,kbmo,kemo,prs_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: prs_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
prs=prs_in
end subroutine set_prs
!------------------
subroutine set_rho(ibmo,iemo,jbmo,jemo,kbmo,kemo,rho_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: rho_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
rho=rho_in
end subroutine set_rho
!------------------
subroutine set_th0(ibmo,iemo,jbmo,jemo,kbmo,kemo,th0_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: th0_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
th0=th0_in
end subroutine set_th0
!------------------
subroutine set_pi0(ibmo,iemo,jbmo,jemo,kbmo,kemo,pi0_in)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(in) :: pi0_in(ibmo:iemo,jbmo:jemo,kbmo:kemo)
pi0=pi0_in
end subroutine set_pi0
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
subroutine get_prs(ibmo,iemo,jbmo,jemo,kbmo,kemo,prs_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: prs_out(ibmo:iemo,jbmo:jemo,kbmo:kemo)
prs_out=prs
end subroutine get_prs
!------------------
subroutine get_rho(ibmo,iemo,jbmo,jemo,kbmo,kemo,rho_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: rho_out(ibmo:iemo,jbmo:jemo,kbmo:kemo)
rho_out=rho
end subroutine get_rho
!------------------
subroutine get_th0(ibmo,iemo,jbmo,jemo,kbmo,kemo,th0_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: th0_out(ibmo:iemo,jbmo:jemo,kbmo:kemo)
th0_out=th0
end subroutine get_th0
!------------------
subroutine get_pi0(ibmo,iemo,jbmo,jemo,kbmo,kemo,pi0_out)
use cm1vars
implicit none
integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo
real, intent(out) :: pi0_out(ibmo:iemo,jbmo:jemo,kbmo:kemo)
pi0_out=pi0
end subroutine get_pi0
!------------------
