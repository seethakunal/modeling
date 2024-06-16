clc;
clear;
close all;
%%
%Givens 
Kt = 1200; 
Vmt = 0.5; 
Vmax = 5; 
Fo = 100; 
Lo = 2; 
Lmt = 2.1;
Lst = 1; 
tau = 0.25;
F_0 = 0.05;
a_0 = 0.05; 
pulsewidth = 0.01;
delay = 1;
t = 0:0.005:5;
stdforce = [];
avgforce = [];
Fstim = [1,2,5,10,15,30,100]; %Hz
for i=1:1:7 
%defining values based on Fstim 
iFstim=Fstim(i);
period = 1/Fstim(i);
if (pulsewidth>period)
    pulsewidth = period;
end 
%ODE solving 
SolverOptions = odeset('MaxStep',pulsewidth/20); %will take longer to run 
[~,AF]=ode23(@CalculateAF,t,[a_0;F_0],SolverOptions,Kt,Vmt,Vmax,Fo,Lo,Lmt,Lst,iFstim,pulsewidth,period,delay);
A = AF(:,1);
F=AF(:,2);
stdforce = [stdforce std(F(500:end))];
avgforce = [avgforce mean(F(500:end))];
figure(1); %plotting f_m
hold on
plot(t,F,'LineWidth',2);
hold off 
figure(2); %plotting A
hold on 
plot(t,A,'LineWidth',2);
hold off 
end 
figure(1); %making the graphs pretty (adding legend and stuff)
hold on
xlabel('Time(s)');
ylabel('F_m (N)');
lgd = legend('1Hz','2Hz','5Hz','10Hz','15Hz','30Hz','100Hz');
lgd.Title.String = 'Fstim=';
hold off
figure(2);
hold on
xlabel('Time(s)');
ylabel('A');
lgd = legend('1Hz','2Hz','5Hz','10Hz','15Hz','30Hz','100Hz');
lgd.Title.String = 'Fstim=';
hold off
figure(3); %plotting ripples
hold on
plot(Fstim, stdforce/max(stdforce))
xlabel('stimulation frequency (hz)')
ylabel('Normalized values')
plot(Fstim, avgforce/max(avgforce))
legend('standard deviation','average force')
hold off 
% ylabel('average muscle force (N)')
%{
1. The source of the ripple is the frequency of the stimulation. 
2. As the frequency as increases, the muscle force also increases. The ripple
exists, but there aren't as many with a low frequency input. Fused tetanus
seems to appear at either 100Hz because there is no dips in force. 
3. Having an activation too low can cause muscle spasms. Therefore, we
would like A to be at least 48 Hz. This is where the standard deviation and
normalized average force meet in the graph. The normalized standard
deviation should be larger or equal to the normalized average force. 
%}
%% 
function [dAF] = CalculateAF(t,AFm,Kt,Vmt,Vmax,Fo,Lo,Lmt,Lst,iFstim,pulsewidth,period,delay)
%givens
tau = 0.25;
%inputs
A = AFm(1);
Fm = AFm(2);
%calc activation dynamics 
% Equation 3
duty = pulsewidth/period*100;
p = square(2*pi*iFstim*(t-delay),duty).*heaviside(t-delay); %The only difference is the first value - w heaviside its .5
%p = square(2*pi*iFstim*(t-delay),duty).*double(t>=delay); %minus one or minus delay?
p(p<0) = 0; %to ensure its never negative 
if A<0
    A=0;
end
%plot(t,p);
%muscle dynamics
Ft = Fm; %because muscle and tendon are in series
LmNorm = (Lmt-Ft/Kt-Lst)/Lo;
LT = 0;
if (LmNorm>=0.42) && (LmNorm <=1.54)
    LT=-3.05*LmNorm^2+5.98*LmNorm-1.96;
end 
VmNorm=Fm/(Fo*A*LT);
%FVinv(VmNorm)
FVinv = 0.995242*exp(13.8817*(VmNorm-1.39214))-0.996815*exp(-3.91442*VmNorm);
if VmNorm<0
    FVinv = -1;
elseif VmNorm>1.4
    FVinv = 1;
end
%derivatives 
dAF(1,1)=(p-A)/tau;
dAF(2,1)=Kt*(Vmt-(Vmax*FVinv));
end 