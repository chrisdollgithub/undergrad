function term_project_423_pt2
clear all
clc;
close all
%Establish Constants
g=9.81;
u0=0;                               
v0=0;
h0=5000;
hmax=5030;    
b=0; %initialize b (ocean floor profile)
%X,Y,t Domain
ni=151;                               
xmax=100000;                      
dx=xmax/(ni-1);
x=0:dx:xmax;
nj=151;                              
ymax=100000;                      
dy=ymax/(nj-1);
y=0:dy:ymax;
%Speed
wavespeed=u0+sqrt(g*(hmax-b));      
%Establish time 
dt = 0.68*dx/wavespeed;
tmax=1500;
t=1:dt:tmax;
c=wavespeed*dt/dx; %courant 
%zero u,v,b,h
u=zeros(length(x),length(y),length(t));v=zeros(length(x),length(y),length(t));b=zeros(length(x),length(y));h=zeros(length(x),length(y),length(t)); 
%Define height (based on Part 1)
h(:,:,1) = h0;                              
h((45000/100000*(length(x)-1)+1):floor(55000/100000*(length(x)-1)+1),...
    (45000/100000*(length(y)-1)+1):floor(55000/100000*(length(y)-1)+1),1) = hmax;
%%%% Part 1)
%Define b (ocean floor profile increases from b=0 until shore starting @ x=20,000)
for i = 1:length(x)
   if x(i)<20000
        b(:,i) = 5000/20000*(20000-x(i)); %%contour of ocean floor within 20,000 m of shore
    elseif x(i)>20001
        b(:,i) = 0;
    end
end
%%%% Part 2)
% Lax 2D
for n=1:(length(t)-1)
    for i=2:(ni-1)
        for j=2:(nj-1)  
            h(i,j,n+1) = ((h(i+1,j,n) + h(i-1,j,n) + h(i,j+1,n) + h(i,j-1,n))/4)- 0.5*(dt/dx)*(u(i,j,n))*((h(i+1,j,n)-b(i+1,j)) - (h(i-1,j,n)-b(i-1,j)))- 0.5*(dt/dy)*(v(i,j,n))*((h(i,j+1,n)-b(i,j+1)) - (h(i,j-1,n)-b(i,j-1))) - 0.5*(dt/dx)*(h(i,j,n)-b(i,j))*(u(i+1,j,n)- u(i-1,j,n)) - 0.5*(dt/dy)*(h(i,j,n)-b(i,j))*(v(i,j+1,n) - v(i,j-1,n));
            u(i,j,n+1) = ((u(i+1,j,n) + u(i-1,j,n) + u(i,j+1,n) + u(i,j-1,n))/4)- 0.5*(dt/dx)*((u(i+1,j,n)^2)/2 - (u(i-1,j,n)^2)/2)- 0.5*(dt/dy)*(v(i,j,n))*(u(i,j+1,n) - u(i,j-1,n)) - 0.5*g*(dt/dx)*(h(i+1,j,n)-h(i-1,j,n));
            v(i,j,n+1) = ((v(i+1,j,n) + v(i-1,j,n) + v(i,j+1,n) + v(i,j-1,n))/4)- 0.5*(dt/dy)*((v(i,j+1,n)^2)/2 - (v(i,j+1,n)^2)/2) - 0.5*(dt/dx)*(u(i,j,n))*(v(i+1,j,n) - v(i-1,j,n)) - 0.5*g*(dt/dy)*(h(i,j+1,n)-h(i,j-1,n));
        end
    end
%BCs (h,x,y)
    h(1,:,n+1) = 2.5*h(2,:,n+1) - 2*h(3,:,n+1) + 0.5*h(4,:,n+1);
    h(length(x),:,n+1) = 2.5*h(ni-1,:,n+1) - 2*h(ni-2,:,n+1) + 0.5*h(ni-3,:,n+1);
    h(:,1,n+1) = 2.5*h(:,2,n+1) - 2*h(:,3,n+1) + 0.5*h(:,4,n+1);
    h(:,length(y),n+1) = 2.5*h(:,nj-1,n+1) - 2*h(:,nj-2,n+1) + 0.5*h(:,nj-3,n+1);
    
    u(1,:,n+1) = 2.5*u(2,:,n+1) - 2*u(3,:,n+1) + 0.5*u(4,:,n+1);
    u(length(x),:,n+1) = 2.5*u(ni-1,:,n+1) - 2*u(ni-2,:,n+1) + 0.5*u(ni-3,:,n+1);
    u(:,1,n+1) = 2.5*u(:,2,n+1) - 2*u(:,3,n+1) + 0.5*u(:,4,n+1);
    u(:,length(y),n+1) = 2.5*u(:,nj-1,n+1) - 2*u(:,nj-2,n+1) + 0.5*u(:,nj-3,n+1);
    %{ 
    v(1,:,n+1)=u(1,:,n+1);
    v(:,1,n+1)=u(:,1,n+1);
    v(length(x),:,n+1)=u(length(x),:,n+1);
    v(:,length(y),n+1)=u(:,length(y),n+1);
    %}
    v(1,:,n+1) = 2.5*v(2,:,n+1) - 2*v(3,:,n+1) + 0.5*v(4,:,n+1);
    v(length(x),:,n+1) = 2.5*v(ni-1,:,n+1) - 2*v(ni-2,:,n+1) + 0.5*v(ni-3,:,n+1);
    v(:,1,n+1) = 2.5*v(:,2,n+1) - 2*v(:,3,n+1) + 0.5*v(:,4,n+1);
    v(:,length(y),n+1) = 2.5*v(:,nj-1,n+1) - 2*v(:,nj-2,n+1) + 0.5*v(:,nj-3,n+1);     
end

%obj = VideoWriter('tsunami.avi');
%obj.Quality = 100;
%obj.FrameRate = 40;
%open(obj);

%Tsunami Propagation (2D Lax's Method)
for index=1:length(t)
    mesh(x,y,h(:,:,index))
    colorbar
    axis ([0 100001 0 100001 4990 5040])
    title ('Term Project Part 2 | AERSP 423')
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Height [m]')
%    f = getframe(gcf)
%    writeVideo(obj,f)

    pause(.0001)
end
%obj.close();
end
