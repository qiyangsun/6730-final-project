% initial conditions
y0 = [1,1];

%parameters
t_end = 100;

options = odeset('RelTol',1e-8,'Maxstep',0.1);
[t,y] = ode15s(@FHNeqns,[0 t_end],y0,options);

%plot numerical solution
figure(1)
subplot(211), plot(t,y(:,1),'b-','LineWidth',2)
hold on
ylabel('x')
set(gca,'Fontsize',14,'LineWidth',1)
subplot(212), plot(t,y(:,2),'b-','LineWidth',2)
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)'), ylabel('y')

%plot numerical solution
figure(2)
plot(y(:,1),y(:,2),'b-','LineWidth',2)
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('x'), ylabel('y')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = FHNeqns(t,y)

a = 0.7; 
b = 0.8;
c = 3.0;

amp = 1;
dur = 0.5;
tstart = 10;
inter = 1.1; 
S=0;
if ((t>tstart && t<tstart+dur) || (t>tstart+inter && t<tstart+dur+inter)) 
    S=amp;
end

dydt = [c*(y(1)-1/3*y(1)^3-y(2)+S); (y(1)+a-b*y(2))/c];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = -0.2; 
b = 0.8;
c = 3.0;

p = [-1/3 ,0, ((b-1)/b), -(a/b)];
result = roots(p);

x = 0.;

e = [c-c*x^2 -c;1/c -b/c];
lambda = eig(e);

