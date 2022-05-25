function a=a_plan(prof,t)
%---input
%prof - motion profile: 'const','trapz','poly'
%t - scalar

%---output
%a - [ax,ay,az]' (3xN matrix)

%load paramters into function workspace. struct called prm
load('Parameters.mat');
T=prm.T; T1=prm.T1; T2=prm.T2; dT=prm.dT;
dX=prm.dX; dY=prm.dY; dZ=prm.dZ;

a=zeros(3,length(t));
for i=1:length(t)
    ti=t(i);
    switch prof
        case 'const'
            a(:,i)=[0,0,0]';
        case 'trapz'
            if ti<T1
                a(:,i)=[dX,dY,dZ]'/(dT*(T-dT));
            elseif ti<T2
                a(:,i)=[0,0,0]';
            else %t<T
                a(:,i)=-[dX,dY,dZ]'/(dT*(T-dT));
            end
        case 'poly'
            a(:,i)=(60*ti/T^3-180*(ti^2)/T^4+120*(ti^3)/T^5)*[dX,dY,dZ]';
    end
end