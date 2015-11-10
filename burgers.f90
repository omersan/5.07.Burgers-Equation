!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Numerical Schemes for solving inviscid Burgers equation
!
!     Lax-Wendroff scheme
!     McCormack scheme
!     Upwind scheme 
!	  WENO scheme
!
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) 
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Nov. 5, 2015
!-----------------------------------------------------------------------------!

program burgers1d
implicit none
integer::i,k,nx,nt,nf,kf
real*8 ::dx,dt,x0,xL,pi,t,Tmax
real*8,allocatable ::u(:),x(:)

!Domain
x0 = 0.0d0 !left
xL = 1.0d0 !right

!number of points
nx = 200

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

!maximum time desired
Tmax = 0.3d0

!number of total points in time
nt = 2000

!number of times to plot
nf = 20

!ploting frequency
kf = nt/nf

!time step
dt = Tmax/dfloat(nt)


!allocate field variable
allocate(u(0:nx))

!initial condition
pi = 4.0d0*datan(1.0d0)
t = 0.0d0
do i=0,nx
u(i) = dsin(2.0d0*pi*x(i))
end do

!boundary conditions: holds for all time
u(0) = 0.0d0
u(nx)= 0.0d0

!Plot initial condition
open(18,file='u.plt')
write(18,*) 'variables ="x","u"'
write(18,100)'zone f=point i=',nx+1,',t="t =',t,'"'
do i=0,nx
write(18,*) x(i),u(i)
end do

!time integration
do k=1,nt
     
	!Lax-Wendroff scheme
    !call LW(nx,dx,dt,u)

    !McCormack scheme
    !call MC(nx,dx,dt,u)
   
    !upwind scheme (flux splitting)
    !call upwind(nx,dx,dt,u)

    !WENO
    call weno(nx,dx,dt,u)
     
    !update t
    t = t+dt 

    !plot the field
    if (mod(k,kf).eq.0) then
	write(18,100)'zone f=point i=',nx+1,',t="t =',t,'"'
	do i=0,nx
	write(18,*) x(i),u(i)
	end do
    end if
    
end do
  

close(18)

100 format(a16,i8,a10,f8.4,a3)
end


!-----------------------------------------------------------------------------!
!Lax-Wendroff scheme
!-----------------------------------------------------------------------------!
subroutine LW(nx,dx,dt,u)
implicit none
integer::nx,i
real*8 ::dx,dt,f1,f0
real*8 ::u(0:nx),up(0:nx)


!predictor step
do i=0,nx-1
	f1 = 0.5d0*u(i+1)*u(i+1)
	f0 = 0.5d0*u(i)*u(i)
	up(i) = 0.5d0*(u(i+1)+u(i)) - 0.5d0*dt/dx*(f1-f0)
end do

!corrector step
do i=1,nx-1
	f1 = 0.5d0*up(i)*up(i)
	f0 = 0.5d0*up(i-1)*up(i-1)
	u(i) = u(i) - dt/dx*(f1-f0)
end do

return
end 

!-----------------------------------------------------------------------------!
!McCormack scheme
!-----------------------------------------------------------------------------!
subroutine MC(nx,dx,dt,u)
implicit none
integer::nx,i
real*8 ::dx,dt,f1,f0
real*8 ::u(0:nx),up(0:nx)

!predictor step
do i=0,nx-1
	f1 = 0.5d0*u(i+1)*u(i+1)
	f0 = 0.5d0*u(i)*u(i)
	up(i) = u(i) - dt/dx*(f1-f0)
end do

!corrector step
do i=1,nx-1
	f1 = 0.5d0*up(i)*up(i)
	f0 = 0.5d0*up(i-1)*up(i-1)
	u(i) = 0.5d0*(u(i)+up(i)) - 0.5d0*dt/dx*(f1-f0)
end do

return
end 


!-----------------------------------------------------------------------------!
!upwind scheme 
!-----------------------------------------------------------------------------!
subroutine upwind(nx,dx,dt,u)
implicit none
integer::nx,i
real*8 ::dx,dt
real*8 ::u(0:nx),v(0:nx),r(1:nx-1)

!uses three stages 3rd order Runge Kutta time integrator
!intermediate variables denoted as v
v(0) = u(0)  !update bc for v
v(nx)= u(nx) !update bc for v

!1st stage
call rhs_upwind(nx,dx,u,r)     
do i=1,nx-1
    v(i) = u(i) + dt*r(i)
end do

!2nd stage	
call rhs_upwind(nx,dx,v,r)
do i=1,nx-1
	v(i) = 0.75d0*u(i) +0.25d0*v(i) + 0.25d0*dt*r(i)
end do

!3rd stage
call rhs_upwind(nx,dx,v,r)
do i=1,nx-1
    u(i) = 1.0d0/3.0d0*u(i) +2.0d0/3.0d0*v(i) + 2.0d0/3.0d0*dt*r(i)
end do

return   
end 



!-----------------------------------------------------------------------------!
!RHS for upwind scheme
!-----------------------------------------------------------------------------!
subroutine rhs_upwind(nx,dx,u,r)
implicit none
integer::nx,i
real*8 ::dx,a,b
real*8 ::u(0:nx),r(1:nx-1)

!first order (using only close to the boundary)
do i=1,nx-1,nx-2    
    a = max(u(i),0.0d0)
    b = min(u(i),0.0d0)
    r(i) = -a*( u(i) - u(i-1))/dx &
           -b*(-u(i) + u(i+1))/dx    
end do

!second order
do i=2,nx-2    
    a = max(u(i),0.0d0)
    b = min(u(i),0.0d0)
    r(i) = -a*( 3.0d0*u(i) - 4.0d0*u(i-1) + u(i-2))/(2.0*dx)  &
           -b*(-3.0d0*u(i) + 4.0d0*u(i+1) - u(i+2))/(2.0*dx)   
end do


return   
end 


!-----------------------------------------------------------------------------!
!WENO scheme 
!-----------------------------------------------------------------------------!
subroutine weno(nx,dx,dt,u)
implicit none
integer::nx,i
real*8 ::dx,dt
real*8 ::u(0:nx),v(0:nx),r(1:nx-1)

!uses three stages 3rd order Runge Kutta time integrator
!intermediate variables denoted as v
v(0) = u(0)  !update bc for v
v(nx)= u(nx) !update bc for v

!1st stage
call rhs_weno(nx,dx,u,r)     
do i=1,nx-1
    v(i) = u(i) + dt*r(i)
end do

!2nd stage	
call rhs_weno(nx,dx,v,r)
do i=1,nx-1
	v(i) = 0.75d0*u(i) +0.25d0*v(i) + 0.25d0*dt*r(i)
end do

!3rd stage
call rhs_weno(nx,dx,v,r)
do i=1,nx-1
    u(i) = 1.0d0/3.0d0*u(i) +2.0d0/3.0d0*v(i) + 2.0d0/3.0d0*dt*r(i)
end do

return   
end 


!-----------------------------------------------------------------------------------!
!Compute right-hand-side (Standard WENO Schemes)
!-----------------------------------------------------------------------------------!
subroutine rhs_weno(nx,dx,u,r)
implicit none
integer:: nx,i
real*8 :: dx
real*8 :: u(0:nx),r(1:nx-1)
real*8 :: q(-1:nx+1)
real*8 :: v1,v2,v3,g

!assing u as q with ghost points 
do i=0,nx
q(i) = u(i)
end do
q(-1)  = 2.0d0*q(0) - q(1)
q(nx+1)= 2.0d0*q(nx) - q(nx-1)


do i=1,nx-1

   	if (u(i).gt.0.0d0) then

     	v1 = (q(i-1) - q(i-2))/dx
     	v2 = (q(i)   - q(i-1))/dx
     	v3 = (q(i+1) - q(i))/dx

  		call w3(v1,v2,v3,g) 
     
     	r(i) = - u(i)*g  
        
   	else
     
     	v1 = (q(i+2) - q(i+1))/dx
     	v2 = (q(i+1) - q(i))/dx
     	v3 = (q(i)   - q(i-1))/dx

   	call w3(v1,v2,v3,g)

     	r(i) = -u(i)*g 
 
   	end if

end do


return
end 

!----------------------------------------------------------------------------------!
!WENO3 
!----------------------------------------------------------------------------------!
subroutine w3(a,b,c,f)
implicit none
real*8 ::a,b,c,f
real*8 ::q1,q2
real*8 ::s1,s2
real*8 ::a1,a2
real*8 ::eps

eps = 1.0d-6

q1 =-0.5d0*a + 1.5d0*b 
q2 = 0.5d0*b + 0.5d0*c

s1 = (b-a)**2
s2 = (c-b)**2

a1 = (1.0d0/3.0d0)/(eps+s1)**2
a2 = (2.0d0/3.0d0)/(eps+s2)**2

f = (a1*q1 + a2*q2)/(a1 + a2)

return
end 









