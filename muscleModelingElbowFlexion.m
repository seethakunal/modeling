tic
clc;clear all;close all;
starting_theta=45;
target_thetas=[0:20:180]; %20
hand_masses=[2.5, 5.5];
ICs=[starting_theta; 0; 0.05; 0.01]; % theta, omega, a, Fm
for j=1:length(hand_masses)
    hand_mass=hand_masses(j);
    
    PA=zeros(length(target_thetas),length(hand_masses));
    for i=1:length(target_thetas)
        target_theta=target_thetas(i);
        
        ub=1;
        lb=0.01;
        temp_amplitude=mean([ub,lb]);
        this_theta=Inf; % was inf
        while (abs(target_theta-this_theta)>1) && ((ub-lb)/(ub+lb)>0.01)
            
            [t,y] = ode23s(@MuscleSystem, [0 10], ICs, [], temp_amplitude, hand_mass);
            this_theta=y(end,1);
            if this_theta > target_theta
                ub=temp_amplitude;
            else
                lb=temp_amplitude;
            end
            temp_amplitude=mean([ub,lb]);
        end
        trajectories{i,j}=[t,y];
        PA(i,j)=temp_amplitude;
    end
    
    figure(j)
    for i=1:length(target_thetas)
        plot(trajectories{i,j}(:,1), trajectories{i,j}(:,2))
        hold on
    end
    legend([num2str(target_thetas')]);
    xlabel('Time (s)')
    ylabel('Thetas (deg)')
    title(['Elbow Angle vs. Time for ' num2str(hand_mass) ' kg'])
    
    figure(3)
    plot(target_thetas, PA(:,j))
    hold on
    xlabel('Theta (deg)')
    ylabel('Threshold Amplitude (normalized)')
    title('Effect of Mass On Threshold')
    legend([num2str(hand_masses') [' kg';' kg']]);
end
toc
%{
1. In the last lab, pulse(t) used the square() function, which was based on stimulus frequency.  In this lab, however, 
pulse(t) is not based on stimulus frequency at all.  Under what condition(s) is this assumption acceptable? 
    
    if stimulus frequency is high
 
2. Why do the 0 and 180 degree lines become flat in Figures 1 and 2? 
    
    boundaries of our system are .1 and 160, and the arm cannot
    extend further beyond that. 
 
3. Why does the threshold amplitude achieve a value of 1.0 at 180 degrees, regardless of the mass of the ball? 
    
    Because, no matter how hard you push, the muscle will never reach that place.  
 
4. What effect does the weight of the ball have on the required amplitude to achieve each angle? 
 
   The PA for the 5.5kg ball is more than then 2.55kg ball, which makes sense
   since more muscle fibers are recruited to lift the ball. 
 
5. With a few exceptions, the arm angle achieved the desired angle only at the very end of the simulation (10 
seconds).  What would you need to do if you wanted the arm to the hold the desired arm angle for an additional 
10 seconds? 
    Keeping the PA where it was at 10 seconds would cause the arm to
    continue to contract. So, to keep the angle the same, the PA would need
    to decrease some amount and stay constant at that level. 
%}
function dy=MuscleSystem(t,y,PA,hand_mass)
theta=y(1);
omega=y(2);
a=y(3); 
Fm=y(4) ;
Io=2.11;
Fo=1800;
Lst=0.1;
Lo=0.25;
Vmax=4.5;
tau=0.25;
Location.forearm=[18,0,0]*.01;
Location.hand=[42,0,0]*.01;
r.forearm=sqrt(sum(Location.forearm.^2)); 
Location.forearm(1)=[r.forearm*sin(deg2rad(theta))];
Location.forearm(2)=-[r.forearm*cos(deg2rad(theta))];
    
r.hand=sqrt(sum(Location.hand.^2)); 
Location.hand(1)=[r.hand*sin(deg2rad(theta))]; 
Location.hand(2)=-[r.hand*cos(deg2rad(theta))];
Force.forearm=[0,1.45,0]*-9.81;
Force.hand=[0,hand_mass,0]*-9.81; 
Location.Origin=[0,30,0]*.01;
Location.Insertion=[8,1,0]*.01;
theta_muscle=rad2deg(atan(Location.Insertion(2)/Location.Insertion(1)));
r.insertion=sqrt(sum(Location.Insertion.^2));
Location.Insertion(1)=[r.insertion*sin(deg2rad(theta+theta_muscle))];
Location.Insertion(2)=-[r.insertion*cos(deg2rad(theta+theta_muscle))];
Velocity=cross(Location.Insertion,[0 0 omega]);
Lmt=sqrt(sum((Location.Insertion-Location.Origin).^2)); 
u=(Location.Insertion-Location.Origin)/Lmt;
Vmt=Velocity*(-u)';
F=Fo*Fm*u;
Moments.forearm=-cross(Force.forearm,Location.forearm); 
Moments.hand=-cross(Force.hand,Location.hand); 
Moments.biceps=-cross(Location.Insertion,F); 
Moments.Total=Moments.forearm+Moments.hand+Moments.biceps; 
dy(1)=y(2);
dy(2)=Moments.Total(3)/Io;
scaling_factor=(1+0.2)*PA/(PA+0.2);
pulse=scaling_factor*heaviside(t-0.5);
dy(3)=(pulse-a)/tau;
dy(4)=MuscleDynamics(t,Fm,Lst,Lo,Lmt,Vmt,Vmax,a);
if theta<0.1
    if omega<0
        dy(1)=0;
    end
    if dy(2)<0
        dy(2)=0;
    end
elseif theta>160
    if omega>0
        dy(1)=0;
    end
    if dy(2)>0
        dy(2)=0;
    end
end
dy=dy';
end
function muscle_dy = MuscleDynamics(t,Fm,Lst,Lo,Lmt,Vmt,Vmax,a)
emax = 0.0325;
Vmtt = Vmt/Lo;
Lmtt = Lmt/Lo;
Lstt = Lst/Lo;
Vmxt = Vmax/Lo;

x1 = (Lmtt-Fm*emax*Lstt-Lstt);
if(x1 >= 0.42 && x1 <= 1.54)
    FL=(-3.0508*x1*x1)+(5.9758*x1)-1.9598;
else
    FL=0;
end
fV=Fm/(a*FL);
if(fV >= 0 && fV <= 1.4)
    Vmn = 0.995242*exp(13.8817*(fV-1.39214))-0.996815*exp(-3.91442*fV);
else
    if(fV<0)
        Vmn=-1;
    else
        Vmn=1;
    end
end
muscle_dy(1)=(Vmtt-Vmxt*Vmn)/(emax*Lstt);
end