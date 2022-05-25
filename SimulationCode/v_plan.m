function v=v_plan(prof,t)
%---input
%prof - motion profile: 'const','trapz','poly'
%t - scalar

%---output
%v - [vx,vy,vz]' (3xN matrix)

%load paramters into function workspace. struct called prm
load('Parameters.mat');
T=prm.T; T1=prm.T1; T2=prm.T2; dT=prm.dT;
dX=prm.dX; dY=prm.dY; dZ=prm.dZ;

v=zeros(3,length(t));
for i=1:length(t)
    ti=t(i);
    switch prof
        case 'const'
            v(:,i)=[dX,dY,dZ]'/T;
        case 'trapz'
            a=[dX,dY,dZ]'/(dT*(T-dT));
            if ti<T1
                v(:,i)=a*ti;
            elseif ti<T2
                v(:,i)=a*T1;
            else %t<T
                v(:,i)=a*T1-a*(ti-T2);
            end
        case 'poly'
            v(:,i)=(30*(ti^2)/T^3-60*(ti^3)/T^4+30*(ti^4)/T^5)*[dX,dY,dZ]';
    end
end