function x=x_plan(prof,t)
%---input
%prof - motion profile: 'const','trapz','poly'
%t - scalar

%---output
%x - [xx,xy,xz]' (3xN matrix)

%load paramters into function workspace. struct called prm
load('Parameters.mat');
T=prm.T; T1=prm.T1; T2=prm.T2; dT=prm.dT;
X0=prm.X(1); Y0=prm.Y(1); Z0=prm.Z(1);
dX=prm.dX; dY=prm.dY; dZ=prm.dZ;

%assume [X0,Y0,Z0]'=[0,0,0]' and add it later
x=zeros(3,length(t));
for i=1:length(t)
ti=t(i);
switch prof
    case 'const'
        v=[dX,dY,dZ]'/T;
        x(:,i)=v*ti;
    case 'trapz'
        a=[dX,dY,dZ]'/(dT*(T-dT));
        if ti<T1
            x(:,i)=0.5*a*(ti^2);
        elseif ti<T2
            x(:,i)=0.5*a*T1^2+a*T1*(ti-T1);
        else %t<T
            x(:,i)=0.5*a*T1^2+a*T1*(T2-T1)+a*T1*(ti-T2)-0.5*a*(ti-T2)^2;
        end
    case 'poly'
        x(:,i)=(10*(ti/T)^3-15*(ti/T)^4+6*(ti/T)^5)*[dX,dY,dZ]';
end
end
x=x+[X0,Y0,Z0]';