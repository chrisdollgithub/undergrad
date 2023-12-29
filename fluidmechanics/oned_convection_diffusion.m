function oned_convection_diffusion
clear all
ni=101;
alpha=1.;
c=10;
x= linspace(0,1,ni);
dx=1./(ni-1)
dt=.00002;
vnm=alpha*dt/dx/dx
cfl=c*dt/dx

%boundary conditions
tn(1)=305.;
tn(ni)=273.;

%initial conditions
  for i = 2:1:ni-1
  tn(i)=273.;
  end

v = VideoWriter ('onedheat_convecton_diffusion.avi');
open (v);  

  for n=1:1:5040
      
% plot current timestep

    if mod(n-1,40) == 0 
    drawnow;
    plot(x,tn,'-m','Color','b','Linewidth',4)
    set(gca,'fontsize', 14);
    xlabel('x')
    ylabel('T(t)')
    legend('t=')
    legend(sprintfc('t= %e ',(n-1)*dt),'Location','northeast')
    title('1D Unsteady Convection-Diffusion')
    grid on
    grid minor
    axis([0 1 270 310])
    M=getframe(gcf)
    writeVideo(v,M);
    hold off
    end
    
    for i = 2:1:ni-1
    tnp1(i)=-cfl*(tn(i)-tn(i-1))+vnm*(tn(i+1)-2.*tn(i)+tn(i-1))+tn(i);
    end
    
%update:
    for i = 2:1:ni-1
    tn(i)=tnp1(i);
    end

  end %end of time loop

close(v);

end