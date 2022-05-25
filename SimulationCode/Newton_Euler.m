function [Rz,tau]=Newton_Euler(prof,t,varargin)
%---input
%prof - motion profile: 'const','trapz','poly'
%t - vector of time values
%varargin input:
%elbows - [t2_elbow,d3_eblow]: accepts -1 or 1 values
%Mchoice - [1 if M=0.5Kg, 2 if M=0.6Kg]

N=length(t);

elbows  = [1,1]; 
Mchoice = 1;
Fe=zeros(6,N);
for ind	= 1:2:length(varargin)
	comm	= lower(varargin{ind});
	switch comm
		case 'elbows'
            elbows	= lower(varargin{ind+1});
        case 'mchoice'
            Mchoice	= lower(varargin{ind+1});
        case 'Fe'
            Fe=lower(varargin{ind+1});
	end
end
load('Parameters.mat');

q=q_plan(prof,t,'elbows',elbows);
q_dot=q_dot_plan(prof,t,'elbows',elbows,'method','jacobian');
q_dot2=q_dot2_plan(prof,t,'elbows',elbows,'method','jacobian');

%% Prepare:
%4 links [world,L1,L2,L3,tool]
u=prm.u; %3x4
r_c=prm.r_c; %3x4
r_e=prm.r_e; %3x4
f_Ritim1=prm.f_Ritim1;
f_Rit0=prm.f_Rit0;
jointType=prm.jointType;
P=cellfun(@(x) strcmp(x,'P'),jointType);
R=cellfun(@(x) strcmp(x,'R'),jointType);

q_dot=[q_dot;zeros(1,N)]; %add fourth link connecting links 4 and 3
q_dot2=[q_dot2;zeros(1,N)]; %add fourth link connecting links 4 and 3

[d_dot,d_dot2]=deal(zeros(4,N));
d_dot(3,:)=q_dot(3,:);
d_dot2(3,:)=q_dot2(3,:);

I=prm.I;
m=prm.m;

g0vec=[0,0,-1]'*prm.g;
%% Loop
[w,alpha,a_c,a_e]=deal(zeros(3,4,N));
[f,M]=deal(zeros(3,4,N));

for s=1:N %run on time samples
    %front recursion
    %First Link: initalize with starting conditions
    r_c(:,3)=[0,0,1]*(q(3,s)-prm.L1/2);
    r_e(:,3)=[0,0,1]*q(3,s);
    
    w(:,1,s)=u(:,1)*q_dot(1,s);
    alpha(:,1,s)=u(:,1)*q_dot2(1,s)+...
        cross(w(:,1,s),u(:,1)*q_dot(1,s));
    a_c(:,1,s)=cross(alpha(:,1,s),r_c(:,1))+...
        cross(w(:,1,s),cross(w(:,1,s),r_c(:,1)));
    a_e(:,1,s)=cross(alpha(:,1,s),r_e(:,1))+...
        cross(w(:,1,s),cross(w(:,1,s),r_e(:,1)));
    for i=2:4 %start looping
        Ritim1=f_Ritim1{i}(q(:,s));
        
        w(:,i,s)=Ritim1'*w(:,i-1,s)+...
            u(:,i)*q_dot(i,s)*R(i);
        alpha(:,i,s)=Ritim1'*alpha(:,i-1,s)+...
            u(:,i)*q_dot2(i,s)*R(i)+...
            cross(w(:,i,s),u(:,i)*q_dot(i,s));
        a_c(:,i,s)=Ritim1'*a_e(:,i-1,s)+...
            cross(alpha(:,i,s),r_c(:,i))+...
            cross(w(:,i,s),cross(w(:,i,s),r_c(:,i)))+...
            u(:,i)*d_dot2(i,s)*P(i)+...
            2*cross(w(:,i,s),u(:,i))*d_dot(i,s)*P(i);
        a_e(:,i,s)=Ritim1'*a_e(:,i-1,s)+...
            cross(alpha(:,i,s),r_e(:,i))+...
            cross(w(:,i,s),cross(w(:,i,s),r_e(:,i)))+...
            u(:,i)*d_dot2(i,s)*P(i)+...
            2*cross(w(:,i,s),u(:,i))*d_dot(i,s)*P(i);
    end
    
    %starting conditions 
    R4t0=f_Rit0{4}(q(:,s));
    f(:,4,s)=m{4}*a_c(:,4,s)+...%force that link 3 applies on link 4
        -R4t0'*Fe(1:3,s)+...%added minus sign here before R4t0
        -m{i}*R4t0'*g0vec;
    M(:,4,s)=cross(r_c(:,4),f(:,4,s))+... %torque that link 3 applies on link 4
        cross(r_e(:,4)-r_c(:,4),-R4t0'*Fe(1:3,s))+...%added minus sign here before R4t0
        I{4}*alpha(:,4,s)+...
        cross(w(:,4,s),I{4}*w(:,4,s));
    
    for i=3:-1:1
       Rip1ti=f_Ritim1{i+1}(q(:,s));
       Rit0=f_Rit0{i}(q(:,s));
       
       f(:,i,s)=m{i}*a_c(:,i,s)+... %force that link im1 applies on i
           Rip1ti*f(:,i+1,s)+...
           -m{i}*Rit0'*g0vec;
       M(:,i,s)=Rip1ti*M(:,i+1,s)+... %moment that link im1 applies on i
           cross(r_c(:,i),f(:,i,s))+...
           cross(r_e(:,i)-r_c(:,i),Rip1ti*f(:,i+1,s))+...
           I{i}*alpha(:,i,s)+...
           cross(w(:,i,s),I{i}*w(:,i,s));
    end
end
%% check tool acceleration
a_e_4_WorldSystem=zeros(3,N);
for s=1:N
    R4t0=f_Rit0{4}(q(:,s));
   a_e_4_WorldSystem(:,s)=R4t0*a_e(:,end,s);
end
a_tool=a_plan(prof,t);
plot((a_e_4_WorldSystem'-a_tool'));
% plot(a_tool')
%% Calculate tau
tau=zeros(3,N);

M_L1=reshape(M(:,1,:),[3,N]);
tau(1,:)=M_L1(3,:);

M_L2=reshape(M(:,2,:),[3,N]);
tau(2,:)=M_L2(1,:);

f_L3=reshape(f(:,3,:),[3,N]);
tau(3,:)=f_L3(3,:);
%% Calculate Rz
f_L0=reshape(f(:,1,:),[3,N]);
Rz=f_L0(3,:);
end