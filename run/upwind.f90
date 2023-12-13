subroutine upwind(u_new,u,a,dx,dt,nx)
implicit none
integer, intent(in) :: nx
real, intent(in) :: dx,dt   
real, intent(in) :: u(nx),a(nx)
real, intent(out) :: u_new(nx)
integer :: i
do i=2,nx-1
    if (a(i) > 0) then
        u_new(i) = u(i) - a(i)*dt/dx*(u(i)-u(i-1))
    else
        u_new(i) = u(i) - a(i)*dt/dx*(u(i+1)-u(i))
    end if
end do

end subroutine upwind



subroutine laxwendroff(u_new,u,a,dx,dt,nx)
    implicit none
    integer :: nx
    real :: dx,dt   
    real :: u(nx),a(nx)
    real, intent(out) :: u_new(nx)
    real :: u2
    integer :: i
    u_new(1)=u(1)
    u_new(nx)=u(nx)
    do i=2,nx-1
        u_new(i)= u(i) - dt/(2.*dx)*(a(i+1)*u(i+1)-a(i-1)*u(i-1)) +&
        dt**2/(2.*dx**2)*((a(i+1)+a(i))/2*(a(i+1)*u(i+1)-a(i)*u(i))-&
        (a(i)+a(i-1))/2*(a(i)*u(i)-a(i-1)*u(i-1)))
    end do
    
end subroutine laxwendroff

subroutine lagrangian_back(xpos_new,u,a,dx,dt,nx,np)
    implicit none
    integer :: nx
    real :: dx,dt   
    real :: u(nx),a(nx)
    real, intent(out) :: xpos_new(nx)
    
    real :: u2
    integer :: i
    u_new(1)=u(1)
    u_new(nx)=u(nx)
    do i=1,nx
        xpos_new(i)=dx/2+(i-1)*dx+u(i)*dt
    end do
    
end subroutine lagrangian_back