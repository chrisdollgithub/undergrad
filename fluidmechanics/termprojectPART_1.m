clear all
clc
close all
%Chris Doll AERSP 423
%%Constants
g=9.81; 
u_0=0;
h0=5000;
hmax=5300;
b=0;
%%X domain stuff
L=100000;
ni=101;
dx = L./(ni-1);
x=1:dx:L;
%%wave speed
wavespeed= u_0 + sqrt(g*(h0-b)); 
%%establish time
dt=.68*dx/wavespeed;
t=1:dt:1500;
lambda=wavespeed*dt/dx;
%initial conditions u
for i=1:1:ni-1
    un(i)=u_0;
    if(i==100)
        un(i)=wavespeed;
    end
end
%initial conditions h
for i = 1:1:ni-1
  hn(i)=h0;
  if(i>=45 && i<=55)
      hn(i)=hmax;
  end
end
v = VideoWriter ('tsunamiSWWE.avi');
open (v); 

for n=1:1:(length(t)-1) 
    
    % plot current timestep
    if mod(n,1) == 0 
    drawnow;
    plot(x,un,'-m','Color','b','Linewidth',3)
    plot(x,hn,'-m','Color','b','Linewidth',3)
    set(gca,'fontsize', 14);
    xlabel('x')
    ylabel('u(t)')
    legend('t=')
    legend(sprintfc('t= %e ',n*dt),'Location','northeast')
    title('Tsunami Model, LAXs Method')
    grid on
    grid minor
    axis([0 100000 4500 5500])
    M=getframe(gcf)
    writeVideo(v,M);
    hold off
    end
    
    %LAX SCHEME    --- indexing error (index exceeds the number of array
    %elements (100)...
     for i=2:1:ni-1
        unp1(i)=un(i+1)+un(i-1)-(un(i)*dt/(2*dx))*(un(i+1)-un(i-1))-(g*dt/(2*dx))*(hn(i+1)-hn(i-1)); 
        hnp1(i)=hn(i+1)+hn(i-1)-(un(i)*dt/(2*dx))*(hn(i+1)-hn(i-1))-(hn(i)*dt/(2*dx))*(un(i+1)-un(i-1));  
     end
  %extrapolation exit bc to avoid reflection:
  unp1(ni)=unp1(ni-1);
  unp1(1)=unp1(2);
  hnp1(ni)=hnp1(ni-1);
  hnp1(1)=hnp1(2);
  %update:
    for i = 1:1:ni
    un(i)=unp1(i);
    end
    for i = 1:1:ni
    hn(i)=hnp1(i);
    end
    %extrapolation exit bc to avoid reflection:
    un(ni-1)=un(ni-2);
    hn(ni-1)=hn(ni-2);
end
    close(v);

























