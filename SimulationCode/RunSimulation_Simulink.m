function [tsim,tau,q_actual,qdot_actual,xyz_actual,q_planned,xyz_planned]=RunSimulation_Simulink()

%tau_time and Fe_wtime are inputs for function workspace
load('Parameters.mat');

hsys=load_system('RobotSimulation');

options=simset('srcworkspace','current',...
    'MaxStep',prm.dt/10,...
    'Solver','ode45');

simOut=sim('RobotSimulation',[],options); %run simulation

tsim=simOut.tout;
tau=simOut.tau.Data';
q_actual=simOut.q_actual.Data';
qdot_actual=simOut.qdot_actual.Data';
xyz_actual=simOut.xyz_actual.Data';
q_planned=simOut.q_planned.Data';
xyz_planned=simOut.xyz_planned.Data';
end