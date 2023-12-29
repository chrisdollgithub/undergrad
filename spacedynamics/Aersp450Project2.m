clear all

%PROJECT 2
%PART A
%PLOTTING TIME HISTORY OF PROJECTION ERROR
w=0.00000184; %rad/s
to=1;
tf=3600;
%%taking time for every 5 seconds from 1 second to 1 hour.
tspan=to:5:tf;
b0=[1;0;0;0];
[t,b] = ode45('myfun',tspan,b0);
%getting qt at each time (true Rodrigues Parameters.
q1=b(:,2)./b(:,1); 
q2=b(:,3)./b(:,1);
q3=b(:,4)./b(:,1);
quat2CRP=q1+q2+q3;
%storing values from loop below
crp2dcm=zeros(3,3,length(i));
c11=zeros(1,length(i));
c12=zeros(1,length(i));
c13=zeros(1,length(i));
c21=zeros(1,length(i));
c22=zeros(1,length(i));
c23=zeros(1,length(i));
c31=zeros(1,length(i));
c32=zeros(1,length(i));
c33=zeros(1,length(i));
%%getting crp2dcm (NEEDED FOR PART C)
for i=1:1:720
q_1=q1(i,:);
q_2=q2(i,:);
q_3=q3(i,:);
C=(1/(1+quat2CRP.'*quat2CRP));
BigQ=[1+(q_1^2)-(q_2^2)-(q_3^2) 2*(q_1*q_2+q_3) 2*(q_1*q_3-q_2); 
      2*(q_2*q_1-q_3) 1-(q_1^2)+q_2^2+q_3^2 2*(q_2*q_3+q_1);
      2*(q_3*q_1-q_2) 2*(q_3*q_2+q_1) 1-(q_1^2)-q_2^2+q_3^2];
%crp2dcm function
crp2dcm(:,:,i)=C*BigQ;
%getting each generated element of the 3x3 dcm
c11(i,:)=crp2dcm(1,1);
c12(i,:)=crp2dcm(1,2);
c13(i,:)=crp2dcm(1,3);
c21(i,:)=crp2dcm(2,1);
c22(i,:)=crp2dcm(2,2);
c23(i,:)=crp2dcm(2,3);
c31(i,:)=crp2dcm(3,1);
c32(i,:)=crp2dcm(3,2);
c33(i,:)=crp2dcm(3,3);
%increment counter
i=i+1;
end
%--------------------------------------------------------------------------
r1=[-0.1517;-0.9669;0.2050];
r2=[-0.8393;0.4494;0.3044];
r3=[-0.0886;-0.5856;0.8000];
r4=[0.8814;-0.0303;0.5202];
f=65; % focal length [mm]
gauss=.0001; %gaussian noise of std 10^-4 mm 
%iterate through the collinearity equations
b1=zeros(3,length(i));
for i=1:1:720
%corrupting the image plane coords with gaussian noise 10^-4mm std
x1(i,:)=-f*((c11(i,:)*r1(1)+c12(i,:)*r1(2)+c13(i,:)*r1(3)))/((c31(i,:)*r1(1)+c32(i,:)*r1(2)+c33(i,:)*r1(3)))+randn*gauss;
x2(i,:)=-f*((c11(i,:)*r2(1)+c12(i,:)*r2(2)+c13(i,:)*r2(3)))/((c31(i,:)*r2(1)+c32(i,:)*r2(2)+c33(i,:)*r2(3)))+randn*gauss;
x3(i,:)=-f*((c11(i,:)*r3(1)+c12(i,:)*r3(2)+c13(i,:)*r3(3)))/((c31(i,:)*r3(1)+c32(i,:)*r3(2)+c33(i,:)*r3(3)))+randn*gauss;
x4(i,:)=-f*((c11(i,:)*r4(1)+c12(i,:)*r4(2)+c13(i,:)*r4(3)))/((c31(i,:)*r4(1)+c32(i,:)*r4(2)+c33(i,:)*r4(3)))+randn*gauss;
y1(i,:)=-f*(c21(i,:)*r1(1)+c22(i,:)*r1(2)+c23(i,:)*r1(3))/(c31(i,:)*r1(1)+c32(i,:)*r1(2)+c33(i,:)*r1(3))+randn*gauss;
y2(i,:)=-f*(c21(i,:)*r2(1)+c22(i,:)*r2(2)+c23(i,:)*r2(3))/(c31(i,:)*r2(1)+c32(i,:)*r2(2)+c33(i,:)*r2(3))+randn*gauss;
y3(i,:)=-f*(c21(i,:)*r3(1)+c22(i,:)*r3(2)+c23(i,:)*r3(3))/(c31(i,:)*r3(1)+c32(i,:)*r3(2)+c33(i,:)*r3(3))+randn*gauss;
y4(i,:)=-f*(c21(i,:)*r4(1)+c22(i,:)*r4(2)+c23(i,:)*r4(3))/(c31(i,:)*r4(1)+c32(i,:)*r4(2)+c33(i,:)*r4(3))+randn*gauss;
%measurement star vector in the body frame
b1(:,i)=1/sqrt((x1(i,:))^2+((y1(i,:))^2+f^2))*[-x1(i,:);-y1(i,:);f];
b2(:,i)=1/sqrt((x2(i,:))^2+((y2(i,:))^2+f^2))*[-x2(i,:);-y2(i,:);f];
b3(:,i)=1/sqrt((x3(i,:))^2+((y3(i,:))^2+f^2))*[-x3(i,:);-y3(i,:);f];
b4(:,i)=1/sqrt((x4(i,:))^2+((y4(i,:))^2+f^2))*[-x4(i,:);-y4(i,:);f];
%time history of projection error for each star
pj1(:,i)=b1(:,i)-crp2dcm(:,:,i)*r1;
pj2(:,i)=b2(:,i)-crp2dcm(:,:,i)*r2;
pj3(:,i)=b3(:,i)-crp2dcm(:,:,i)*r3;
pj4(:,i)=b4(:,i)-crp2dcm(:,:,i)*r4;
%magnitudes (projection error)
norm1(i,:)=norm(pj1(:,i));
norm2(i,:)=norm(pj2(:,i));
norm3(i,:)=norm(pj3(:,i));
norm4(i,:)=norm(pj4(:,i));
%increment counter
i=i+1;
end
%plot time history of projection error
figure(1)
plot(norm1) 
hold on
plot(norm2) 
hold on
plot(norm3) 
hold on
plot(norm4) 
hold on
title('Time History of Projection Error')
xlabel('Time (s)')
ylabel('Projection Error')
%--------------------------------------------------------------------------
%PART B
%OLAE ALGORITHM
qest1=zeros(3,length(i));
qest2=zeros(3,length(i));
qest3=zeros(3,length(i));
qest4=zeros(3,length(i));
for i=1:1:720
%getting d and s vectors, and vertcat-ing them in order to get column vec 
d1(:,i)=b1(:,i)-r1;
d1stack=vertcat(d1(:));
d2(:,i)=b2(:,i)-r2;
d2stack=vertcat(d2(:));
d3(:,i)=b3(:,i)-r3;
d3stack=vertcat(d3(:));
d4(:,i)=b4(:,i)-r4;
d4stack=vertcat(d4(:));
%getting s vector(s) and then getting the s skew matrices for each; stack
s1(:,i)=r1+b1(:,i);
s2(:,i)=r2+b2(:,i);  
s3(:,i)=r3+b3(:,i);  
s4(:,i)=r4+b4(:,i);
stil1(:,:,i)=[0 -s1(3,i) s1(2,i); s1(3,i) 0 -s1(1,i); -s1(2,i) s1(1,i) 0];
stil1_=num2cell(stil1,[1 2]);
sstack1=vertcat(stil1_{:});
stil2(:,:,i)=[0 -s2(3,i) s2(2,i); s2(3,i) 0 -s2(1,i); -s2(2,i) s2(1,i) 0];
stil2_=num2cell(stil2,[1 2]);
sstack2=vertcat(stil2_{:});
stil3(:,:,i)=[0 -s3(3,i) s3(2,i); s3(3,i) 0 -s3(1,i); -s3(2,i) s3(1,i) 0];
stil3_=num2cell(stil3,[1 2]);
sstack3=vertcat(stil3_{:});
stil4(:,:,i)=[0 -s4(3,i) s4(2,i); s4(3,i) 0 -s4(1,i); -s4(2,i) s4(1,i) 0];
stil4_=num2cell(stil4,[1 2]);
sstack4=vertcat(stil4_{:});
%using OLAE formula to get estimated Rodrigues Parameters
%I'm assuming for each star....
qest1(:,i)=pinv(sstack1)*d1stack;
qest2(:,i)=pinv(sstack2)*d2stack;
qest3(:,i)=pinv(sstack3)*d3stack;
qest4(:,i)=pinv(sstack4)*d4stack;
%increment counter
i=i+1;
end
%%plotted each component of each q (est and true) vector vs time
figure(2)
plot(qest1(1,:))
hold on
plot(qest1(2,:))
hold on
plot(qest1(3,:))
hold on
plot(quat2CRP(1,:))
hold on
legend('qest1(1,:)','qest1(2,:)','qest1(3,:)','quat2crp')
title('Star 1: Estimated RP vs True RP(quat2crp)')
xlabel('Time (s)')
ylabel('RPs')
figure(3)
plot(qest2(1,:))
hold on
plot(qest2(2,:))
hold on
plot(qest2(3,:))
hold on
plot(quat2CRP(1,:))
hold on
legend('qest2(1,:)','qest2(2,:)','qest2(3,:)','quat2crp')
title('Star 2: Estimated RP vs True RP(quat2crp)')
xlabel('Time (s)')
ylabel('RPs')
figure(4)
plot(qest3(1,:))
hold on
plot(qest3(2,:))
hold on
plot(qest3(3,:))
hold on
plot(quat2CRP(1,:))
hold on
legend('qest3(1,:)','qest3(2,:)','qest3(3,:)','quat2crp')
title('Star 3: Estimated RP vs True RP(quat2crp)')
xlabel('Time (s)')
ylabel('RPs')
figure(5)
plot(qest4(1,:))
hold on
plot(qest4(2,:))
hold on
plot(qest4(3,:))
hold on
plot(quat2CRP(1,:))
hold on
legend('qest4(1,:)','qest4(2,:)','qest4(3,:)','quat2crp')
title('Star 4: Estimated RP vs True RP(quat2crp)')
xlabel('Time (s)')
ylabel('RPs')
%--------------------------------------------------------------------------
%PART C
%using estimated rodrigues paramters calculated from part b...
%calculate the estimated, error DCM and write DCM2CRP
q11=zeros(1,length(i));
for i=1:1:720
%estimated DCM Cbn(estdcm)
%Star 1
const1(i,:)=1/((1+qest1(:,i)).'*qest1(:,i));
BigQ1(:,:,i)=[1+(qest1(1,i).^2)-(qest1(2,i).^2)-(qest1(3,i).^2) 2*(qest1(1,i)*qest1(2,i)+qest1(3,i)) 2*(qest1(1,i)*qest1(3,i)-qest1(2,i)); 
       2*(qest1(2,i)*qest1(1,i)-qest1(3,i)) 1-(qest1(2,i).^2)+qest1(2,i).^2+qest1(3,i).^2 2*(qest1(2,i)*qest1(3,i)-qest1(1,i));
       2*(qest1(3,i)*qest1(1,i)-qest1(2,i)) 2*(qest1(3,i)*qest1(2,i)+qest1(1,i)) 1-qest1(1,i).^2-qest1(2,i).^2+qest1(3,i).^2];
estdcm1(:,:,i)=const1(i,:)*BigQ1(:,:,i);
%Star 2
const2(i,:)=1/((1+qest2(:,i)).'*qest2(:,i));
BigQ2(:,:,i)=[1+(qest2(1,i).^2)-(qest2(2,i).^2)-(qest2(3,i).^2) 2*(qest2(1,i)*qest2(2,i)+qest2(3,i)) 2*(qest2(1,i)*qest2(3,i)-qest2(2,i)); 
       2*(qest2(2,i)*qest2(1,i)-qest2(3,i)) 1-(qest2(2,i).^2)+qest2(2,i).^2+qest2(3,i).^2 2*(qest2(2,i)*qest2(3,i)-qest2(1,i));
       2*(qest2(3,i)*qest2(1,i)-qest2(2,i)) 2*(qest2(3,i)*qest2(2,i)+qest2(1,i)) 1-qest2(1,i).^2-qest2(2,i).^2+qest2(3,i).^2];
estdcm2(:,:,i)=const2(i,:)*BigQ2(:,:,i);
%Star 3
const3(i,:)=1/((1+qest3(:,i)).'*qest3(:,i));
BigQ3(:,:,i)=[1+(qest3(1,i).^2)-(qest3(2,i).^2)-(qest3(3,i).^2) 2*(qest3(1,i)*qest3(2,i)+qest3(3,i)) 2*(qest3(1,i)*qest3(3,i)-qest3(2,i)); 
       2*(qest3(2,i)*qest3(1,i)-qest3(3,i)) 1-(qest3(2,i).^2)+qest3(2,i).^2+qest3(3,i).^2 2*(qest3(2,i)*qest2(3,i)-qest3(1,i));
       2*(qest3(3,i)*qest3(1,i)-qest3(2,i)) 2*(qest3(3,i)*qest3(2,i)+qest3(1,i)) 1-qest3(1,i).^2-qest3(2,i).^2+qest3(3,i).^2];
estdcm3(:,:,i)=const3(i,:)*BigQ3(:,:,i);
%Star 4
const4(i,:)=1/((1+qest4(:,i)).'*qest4(:,i));
BigQ4(:,:,i)=[1+(qest4(1,i).^2)-(qest4(2,i).^2)-(qest4(3,i).^2) 2*(qest4(1,i)*qest4(2,i)+qest4(3,i)) 2*(qest4(1,i)*qest4(3,i)-qest4(2,i)); 
       2*(qest4(2,i)*qest4(1,i)-qest4(3,i)) 1-(qest4(2,i).^2)+qest4(2,i).^2+qest4(3,i).^2 2*(qest4(2,i)*qest4(3,i)-qest4(1,i));
       2*(qest4(3,i)*qest4(1,i)-qest4(2,i)) 2*(qest4(3,i)*qest4(2,i)+qest4(1,i)) 1-qest4(1,i).^2-qest4(2,i).^2+qest4(3,i).^2];
estdcm4(:,:,i)=const4(i,:)*BigQ4(:,:,i);
%Getting the error direction cosine matrix for each star
delC1(:,:,i)=estdcm1(:,:,i)*((estdcm1(:,:,i).'));
delC2(:,:,i)=estdcm2(:,:,i)*((estdcm2(:,:,i).'));
delC3(:,:,i)=estdcm3(:,:,i)*((estdcm3(:,:,i).'));
delC4(:,:,i)=estdcm4(:,:,i)*((estdcm4(:,:,i).'));
%Getting error Rodrigues Parameters corresponding to error dcm(skew matrix)
%using the inverse transformation of the Cayley transform...
dcm2crp1(:,:,i)=inv(eye(3)+delC1(:,:,i))*(eye(3)-delC1(:,:,i));
dcm2crp2(:,:,i)=inv(eye(3)+delC2(:,:,i))*(eye(3)-delC2(:,:,i));
dcm2crp3(:,:,i)=inv(eye(3)+delC3(:,:,i))*(eye(3)-delC3(:,:,i));
dcm2crp4(:,:,i)=inv(eye(3)+delC4(:,:,i))*(eye(3)-delC4(:,:,i));
%getting RP from skew matrix [Q] for each star...
q11(:,i)=dcm2crp1(3,2,i);
q12(:,i)=dcm2crp1(1,3,i);
q13(:,i)=dcm2crp1(2,1,i);
q21(:,i)=dcm2crp2(3,2,i);
q22(:,i)=dcm2crp2(1,3,i);
q23(:,i)=dcm2crp2(2,1,i);
q31(:,i)=dcm2crp3(3,2,i);
q32(:,i)=dcm2crp3(1,3,i);
q33(:,i)=dcm2crp3(2,1,i);
q41(:,i)=dcm2crp4(3,2,i);
q42(:,i)=dcm2crp4(1,3,i);
q43(:,i)=dcm2crp4(2,1,i);
%increment counter
i=i+1;
end
%plot time history of error Rodrigues parameters
figure(6)
plot(q11(1,:))
hold on
plot(q12(1,:))
hold on
plot(q13(1,:))
hold on
legend('q11','q12','q13')
title('Star 1: Time History of Error Rodrigues Parameters')
xlabel('Time (s)')
ylabel('error RPs')
figure(7)
plot(q21(1,:))
hold on
plot(q22(1,:))
hold on
plot(q23(1,:))
hold on
legend('q21','q22','q23')
title('Star 2: Time History of Error Rodrigues Parameters')
xlabel('Time (s)')
ylabel('error RPs')
figure(8)
plot(q31(1,:))
hold on
plot(q32(1,:))
hold on
plot(q33(1,:))
hold on
legend('q31','q32','q33')
title('Star 3: Time History of Error Rodrigues Parameters')
xlabel('Time (s)')
ylabel('error RPs')
figure(8)
plot(q41(1,:))
hold on
plot(q42(1,:))
hold on
plot(q43(1,:))
hold on
legend('q41','q42','q43')
title('Star 4: Time History of Error Rodrigues Parameters')
xlabel('Time (s)')
ylabel('error RPs')
%%Some of the estimate Rodrigues parameters seemed to map the same behavior
%%as the true RP relatively accurately.  Some of them were off by a more
%%noticable margin but behavior-wise the estimated RP acted similar to the
%%true RP. 
%--------------------------------------------------------------------------
%PART D
%Using Davenport's q method to compute estimated quaternion parameter at
%each time.  Plot time histories of b est vs b true. 
%I'm going to make wk=1 because I don't know what it should be (non-negative)
BB1=zeros(3,3,length(i));
for i=1:1:720
%Computing the B sum for the k matrix
BB1(:,:,i)=(b1(:,i))*(r1.');
BigB1=sum(BB1,3);
BB2(:,:,i)=(b2(:,i))*(r2.');
BigB2=sum(BB2,3);
BB3(:,:,i)=(b3(:,i))*(r3.');
BigB3=sum(BB3,3);
BB4(:,:,i)=(b4(:,i))*(r4.');
BigB4=sum(BB4,3);
%Computing sigma
sigma1=trace(BigB1);sigma2=trace(BigB2);sigma3=trace(BigB3);sigma4=trace(BigB4);
%Computing S
SS1=BigB1+BigB1.';SS2=BigB2+BigB2.';SS3=BigB3+BigB3.';SS4=BigB4+BigB4.';
%Computing Z
Z1=[BigB1(2,3)-BigB1(3,2) BigB1(3,1)-BigB1(1,3) BigB1(1,2)-BigB1(2,1)].';
Z2=[BigB2(2,3)-BigB2(3,2) BigB2(3,1)-BigB2(1,3) BigB2(1,2)-BigB2(2,1)].';
Z3=[BigB3(2,3)-BigB3(3,2) BigB3(3,1)-BigB3(1,3) BigB3(1,2)-BigB3(2,1)].';
Z4=[BigB4(2,3)-BigB4(3,2) BigB4(3,1)-BigB4(1,3) BigB4(1,2)-BigB4(2,1)].';
%Finally calculate k matrix
k1=[sigma1 Z1.';Z1 (SS1-sigma1*eye(3))];
k2=[sigma2 Z2.';Z2 (SS2-sigma2*eye(3))];
k3=[sigma3 Z3.';Z3 (SS3-sigma3*eye(3))];
k4=[sigma4 Z4.';Z4 (SS4-sigma4*eye(3))];
%get the eigenvalues and eigenvectors
%normalize eigenvectors corresponding to maximum eigenvalue
[b1eVECT,b1eVAL]=eig(k1);
b1eVECT_=b1eVECT(:,3);
b1eVECT=b1eVECT_/(sqrt(b1eVECT_(1,:)^2+b1eVECT_(2,:)^2+b1eVECT_(3,:)^2+b1eVECT_(4,:)^2));
[b2eVECT,b2eVAL]=eig(k2);
b2eVECT_=b2eVECT(:,3);
b2eVECT=b2eVECT_/(sqrt(b2eVECT_(1,:)^2+b2eVECT_(2,:)^2+b2eVECT_(3,:)^2+b2eVECT_(4,:)^2));
[b3eVECT,b3eVAL]=eig(k3);
b3eVECT_=b3eVECT(:,3);
b3eVECT=b3eVECT_/(sqrt(b3eVECT_(1,:)^2+b3eVECT_(2,:)^2+b3eVECT_(3,:)^2+b3eVECT_(4,:)^2));
[b4eVECT,b4eVAL]=eig(k4);
b4eVECT_=b4eVECT(:,3);
b4eVECT=b4eVECT_/(sqrt(b4eVECT_(1,:)^2+b4eVECT_(2,:)^2+b4eVECT_(3,:)^2+b4eVECT_(4,:)^2));
%b(i)eVECT are the estimated quaternion parameters
%increment counter
i=i+1;
end




