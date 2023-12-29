clear all
%Part A
%Calculating the transfer orbit semimajor axis a_min that results in a
%minimum energy transfer orbit. 
mu=1; %AU^3/TU^2
delf=70*pi/180; %convert from degrees to rads
r_earth=149.5*10^6; %km
r_mars=227.8*10^6; %km
r1=r_earth/r_earth; %AU ...normalize r1 and r2 in order to make life easier
r2=r_mars/r_earth; %AU
c=sqrt(r1^2+r2^2-2*r1*r2*cos(delf)); %AU ...calculating chord length
a_min=1/4*(r1+r2+c) %AU ...finally getting min energy xfer semimajor axis
%Part b -------------------------------------------------------------------
%To calculate semi-latus rectum p and eccentricity e,
%I'll use the k,l, and m expressions to calculate p_min, and use that to
%get e_min
k=r1*r2*(1-cos(delf));
l=r1+r2;
m=r1*r2*(1+cos(delf));
num=2*a_min*k*l-k*m+k*sqrt(m*(8*a_min^2-4*a_min*l+m));
den=2*(a_min*l^2-2*a_min*m);
p_min=num/den;
p_min
e_min=sqrt(1-p_min/a_min);
e_min
%Part c -------------------------------------------------------------------
%In order to calculate the initial velocity vector v1 of the xfer orbit,
%I'll use the values calculated in part b.
%F and G solutions using perifocal coords for r1 and r2 vectors...
F=1-r2/p_min*(1-cos(delf));
G=(r1*r2/sqrt((mu*p_min))*sin(delf));
Gdot=1-r1/p_min*(1-cos(delf));
r1vec=[r1 0];
r2vec=[r2*cos(delf) r2*sin(delf)]; 
r1dot=(1/G)*[-F*r1vec+r2vec]; %in terms of p and q unit vectors in that order
r2dot=(1/G)*[-r1vec+Gdot*r2vec];
v1=r1dot %AU/TU
v2=r2dot; %AU/TU
%Part d -------------------------------------------------------------------
%Determine the time of flight delT.
%To begin, use the Kepler Equation to solve for f1 and f2 (using 
%change of true anomaly to solve)...
f1=acos(1/e_min*(p_min/r1-1)); %rad
f2=acos(1/e_min*(p_min/r2-1)); %rad
%Sanity check: 
f2_check=2*pi-acos(1/e_min*(p_min/r2-1)); %this value comes out to be 3.405 rad
check=(f2_check-f1)*180/pi; %...comes out to be 70 degrees, which is our delf.
%Therefore our original f2 is incorrect(f2-f1 =/= 70 deg) so we will use 
%the value of f2_check as our f2...
f2=f2_check; %rad
%using the found values of f1 and f2, go find corresponding 
%eccentric anomaly using the transformation between f and E...
E1=2*atan(tan(f1/2)/(sqrt((1+e_min)/(1-e_min)))); %rad
E2=2*atan(tan(f2/2)/(sqrt((1+e_min)/(1-e_min)))); %rad
%because of the nature of the arctan function in Matlab, atan will return
%the angle between -pi to pi, however we are looking for the angle between
%zero and 2pi so we are working in the same quadrant.
if E2 < 0
    E2 = E2+2*pi;
end
%then calculate mean anomaly
M1=E1-e_min*sin(E1);
M2=E2-e_min*sin(E2);
%Now that we have values for both mean anomalies, we can use the 
%(M2-M1)*sqrt(a_min^3) relation to get our time of flight delT...
delT=(M2-M1)*sqrt(a_min^3) %TU (canonical time units) where mu=AU^3/TU^2
%Part E ------------------------------------------------------------------
%To begin, vary the semimajor axis of xfer orbit from a_min to a_min+5AU
%Then we will run through Parts B-D again with our new a's...
%First, let's generate values of a from a_min to a_min+5 
a=[];
for i=0:0.01:5
 a=[a a_min+i];
end
%Re-do Part B-D
for i=1:1:501
%Re-do Part B
num(:,i)=2*a(:,i)*k*l-k*m+k*sqrt(m*(8*a(:,i)^2-4*a(:,i)*l+m));
den(:,i)=2*(a(:,i)*l^2-2*a(:,i)*m);
p_min(:,i)=num(:,i)/den(:,i);
e_min(:,i)=sqrt(1-p_min(:,i)/a(:,i));
%Re-do Part C 
F(:,i)=1-r2/p_min(:,i)*(1-cos(delf));
G(:,i)=(r1*r2/sqrt((mu*p_min(:,i)))*sin(delf));
Gdot(:,i)=1-r1/p_min(:,i)*(1-cos(delf));
r1vec=[r1 0];
r2vec=[r2*cos(delf) r2*sin(delf)]; 
r1dot(i,:)=(1/G(:,i))*[-F(:,i)*r1vec+r2vec]; %in terms of p and q unit vectors in that order
r2dot(i,:)=(1/G(:,i))*[-r1vec+Gdot(:,i)*r2vec];
%finally getting v1 and v2 at all data points
v1(i,:)=r1dot(i,:);%AU/TU
v2(i,:)=r2dot(i,:);  %AU/TU
%Re-do Part D
%change of true anomaly to solve...
f1(:,i)=acos(1/e_min(:,i)*(p_min(:,i)/r1-1)); %rad
f2(:,i)=acos(1/e_min(:,i)*(p_min(:,i)/r2-1)); %rad
%using the found values of f1 and f2, go find corresponding 
%eccentric anomaly using the transformation between f and E...
E1(:,i)=2*atan(tan(f1(:,i)/2)/(sqrt((1+e_min(:,i))/(1-e_min(:,i))))); %rad
E2(:,i)=2*atan(tan(f2(:,i)/2)/(sqrt((1+e_min(:,i))/(1-e_min(:,i))))); %rad
%because of the nature of the arctan function in Matlab, atan will return
%the angle between -pi to pi, however we are looking for the angle between
%zero and 2pi so we are working in the same quadrant.
if E1 < 0
    E1=E1+2*pi;
end
if E2 < 0
    E2=E2+2*pi;
end
%then calculate mean anomaly
M1(:,i)=E1(:,i)-e_min(:,i)*sin(E1(:,i));
M2(:,i)=E2(:,i)-e_min(:,i)*sin(E2(:,i));
%Now that we have values for both mean anomalies, we can use the 
%(M2-M1)*sqrt(a_min^3) relation to get our time of flight delT...
delT(:,i)=(M2(:,i)-M1(:,i))*sqrt((a(:,i)^3)/mu); %TU (canonical time units) where mu=AU^3/TU^2
%Now to calculate deltaV's...
%Use previous Lambert solution to get v1_sc (spacecraft velocity wrt sun at
%departing point)
r0=6371+600; %[km] earths radius + parking orbit altitude
v1_sc(i,:)=norm(v1(i,:)); %r1dot from lambert solution is s/c's velocity wrt sun 
%at departing point (magnitude)
v1_sc(i,:)=v1_sc(i,:)*(149.5*10^6)/(5.0177*10^6);%convert from AU/TU to km/s
%circular earth park5kming orbit velocity:
mu_sun=1.326*10^11; %km^3/s^2 (grav const of sun)
v_c=sqrt(mu_sun/r_earth);
%required earth relative velocity the s/c requires as it approaches SOI...
v_1(i,:)=v1_sc(i,:)+v_c; 
%because the trajectory will have asymptotically approached hyperbolic 
%asymptote, we make r in the vis-viva equation set to infinity and solve
%for semimajor axis a at the generated values of v_1:
mu_earth=398600.4349; %km^3/s^2
a_hb(i,:)=-1*(mu_earth)/((v_1(i,:))^2);
%Obtain earth relative speed v0:
v0(i,:)=sqrt((2*mu_earth/r0)-(mu_earth/a_hb(i,:)));
%deltaV (burn) to enter the departure hyperbolic orbit is calculated as:
delV(i,:)=(v0(i,:))-v_c;
end
%Part E--------------------------------------------------------------------
figure(1)
plot(a,delT)
title('Plot of T.O.F. as a Function of Semimajor Axis')
xlabel('Semimajor Axis a (AU)')
ylabel('time of flight delT (TU)')
figure(2)
plot(a,delV)
title('Plot of Delta V as a Function of Semimajor Axis')
xlabel('Semimajor Axis a (AU)')
ylabel('Delta V (km/s)')
%COMMENTS: It appears that the plot of tof vs a has 1/x type behavior, so
%time of flight decreases as the semimajor axis of the orbit increases
%until tof approaches infinity asymptotically.  The plot of delta V appears
%to show behavior where it slopes off asymptotically as semi major axis
%increases.
%Part F--------------------------------------------------------------------
%Pick a transfer orbit which corresponds to a transfer time of 3 months...
%Number of seconds in 3 months is approximately 3*2.628*10^6 s 
% => 7.884*10^6 seconds
%Get value of TU from 3 months => 
three_months_TU=(7.884*10^6)/(5.0177*10^6);
%the above value comes out to be roughly 1.57, so look at the value of plot
%1 @ TU=1.57, which appears to be a=1.6184.
figure(3)
x=0;
y=0;
xmin=1.2619;
theta=0:pi/50:2*pi;
hold on
%plot earth
xearth=r1*cos(theta);
yearth=r1*sin(theta);
plotearth=plot(xearth,yearth);
%plot mars
xmars=r2*cos(theta);
ymars=r2*sin(theta);
plotmars=plot(xmars,ymars);
%plot minimum energy transfer orbit
aplot=1.0084;
bplot=0.81834;
xmin=xmin+aplot*cos(theta);
ymin=y+bplot*sin(theta);
plotmin=plot(xmin,ymin);
%plot transfer orbit with time of 3 months
axplot=1.6184;
bxplot=1.4691;
xxfer=xmin+axplot*cos(theta);
yxfer=y+bxplot*sin(theta);
plotmin=plot(xxfer,yxfer);
%Couldn't figure out how to do the plots correctly..



















