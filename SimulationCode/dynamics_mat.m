function [H,C,G]=dynamics_mat(q,qdot,varargin)
%---input
%q - 3xN [t1,t2,d3]' 
%qdot - 3xN [t1dot,t2dot,d3dot]' 

%varargin input:
%Mchoice - [1 if M=0.5Kg, 2 if M=0.6Kg]

%---output
%H,C,G - [3x3xN] matrices of dynamic equations where third dimension
%indicates progress in time
Mchoice=1;
for ind	= 1:2:length(varargin)
	comm	= lower(varargin{ind});
	switch comm
		case 'mchoice'
			Mchoice	= lower(varargin{ind+1});
	end
end

%load parameters
load('Parameters.mat');

[H,C]=deal(zeros(3,3,size(q,2)));
G=zeros(3,size(q,2));
if Mchoice==1
    for i=1:size(q,2)
        H(:,:,i)=prm.f_HdynM1(q(:,i));
        C(:,:,i)=prm.f_CdynM1([q(:,i);qdot(:,i)]);
        G(:,i)=prm.f_GdynM1(q(:,i));
    end
else %Mchoice==2
    for i=1:size(q,2)
        H(:,:,i)=prm.f_HdynM2(q(:,i));
        C(:,:,i)=prm.f_CdynM2([q(:,i);qdot(:,i)]);
        G(:,i)=prm.f_GdynM2(q(:,i));
    end
end