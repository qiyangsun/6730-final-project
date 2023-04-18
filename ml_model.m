%file Morris_Lecar_cable.m
%Model of AP propagation

global L;
global t_end;
global D;
global x_on;
global V_on;

L = 30; %length of neuron (in cm)
t_end = 300; %time to simulate over (in ms)
D = 0.0375; %diffusion coeffcient cm^2/ms
x_on = 1; %cm, length over which V is elevated initially
V_on = 20; %mV; value V is elevated to initially

%set up solution mesh
m = 0;
Nx = 200;
Nt = 200;
x = linspace(0,L,Nx);
t = linspace(0,t_end,Nt);

%solve the partial differential equations
sol = pdepe(m,@Morris_Lecar_cable_eqn,@Morris_Lecar_cable_initial,@Morris_Lecar_cable_bc,x,t);

V = sol(:,:,1);
w = sol(:,:,2);

%plot results
figure(1)
surf(x,t,V, 'EdgeColor','none')    
xlabel('Position (cm)')
ylabel('Time (ms)')
zlabel('Membrane Voltage (mV)')
shading interp
set(gca,'fontsize',14)
figure(2)
surf(x,t,V, 'EdgeColor','none')    
xlabel('Position (cm)')
ylabel('Time (ms)')
zlabel('Membrane Voltage (mV)')
shading interp
set(gca,'fontsize',14)
view(2)
hold on


%plot spatial profiles at a few different time points
% figure(2)
% set(gca,'fontsize',14)
% hold on
% plot(x, V(1,:),'k', 'linewidth', 2);
% plot(x, V(50,:),'k--', 'linewidth', 2);
% plot(x, V(150,:),'k:', 'linewidth', 2);
% %axis([0 15 -80 40])
% xlabel('Position (cm)')
% ylabel('Membrane voltage (mV)')
% legend('t=0 ms', 't=50 ms', 't=150 ms')
% 
% figure(2)
% plot(t,V(:,50),'k')
% hold on
% plot(t,V(:,100),'r')
% hold on
% plot(t,V(:,120),'b')
% xlabel('Time (ms)')
% ylabel('Membrane Voltage (mV)')
% set(gca,'fontsize',14)

%set equations
function [c,f,s] = Morris_Lecar_cable_eqn(x,t,u,DuDx)

global D;

Cap=20;
V_Nernst_K=-84;
gbar_K=8;
V_Nernst_Ca=120;
gbar_Ca=4.4;
gbar_leak=0.5;
V_Nernst_leak=-60;
v1=-1.2;
v2=18;
v3=2;
v4=30;

m_inf = @(x) ( 0.5*(1+tanh((x-v1)/v2)) );
w_inf = @(x) ( 0.5*(1+tanh((x-v3)/v4)) );
tau_w=0.8/.04;

c = [1;1];
f = [D;0].*DuDx;
s = [(1/Cap)*(-gbar_Ca*m_inf(u(1))*(u(1)-V_Nernst_Ca)-gbar_K*u(2)*(u(1)-V_Nernst_K)...
        -gbar_leak*(u(1)-V_Nernst_leak)); 
    (w_inf(u(1))-u(2))/tau_w];

end


%set initial condition
function value = Morris_Lecar_cable_initial(x)

global x_on;
global V_on;
% 
if x<=x_on || x>=30-x_on
    value = [V_on; 0.004];
else
    value = [-61.9; 0.004];
end
%value = [V_on; 0.004];
end


%set boundary conditions
function [pl,ql,pr,qr] = Morris_Lecar_cable_bc(xl,ul,xr,ur,t)

pl = [0;0];
ql = [1;1];
pr = [0;0];
qr = [1;1];

end
