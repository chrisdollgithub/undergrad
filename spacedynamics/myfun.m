%%function for project 2 angl velocity 
function bdot=myfun(t,b)
w=0.00000184; %rad/s
w1=w*sin(w*t);
w2=w*cos(w*t);
w3=w;
BigW =[0, -w1, -w2, -w3;
       w1, 0, w3, -w2;
       w2, -w3, 0, w1;
       w3, w2, -w1, 0];
bdot=0.5*BigW*b;
end

