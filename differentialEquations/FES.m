function [dy] = FES(t, y, m)
% map:
% y = y(1)
% dy(1) = y' = y(2)
% dy(2) = y''

% Force external
Fext = 80*heaviside(t-8)*exp(-(t-8)/2);

% B corresponds to data for biceps, T is for triceps, m is mass
ActB = zeros(size(t));
ActB(t > 10) = 0.35;
ActT = .15;

% resting length
Lrest = 18*(ActB-ActT);

% spring constants
kT = ActT*40;
kB = ActB*50;

% damping constants
cT = 1 + kT/6;
cB = 1 + kB/7;

% calculates the differential equations in the output dy
dy(1) = y(2);
dy(2) = -((cB+cT)/m)*y(2) - ((kB+kT)/m)*y(1) + (kT*Lrest + kB*Lrest)/m + Fext/m;

dy = dy'; % converts dy to column vector by taking transpose

% solver times 
if isempty(evalin('base','SolverTimes' ))
    evalin('base',['SolverTimes(1) = ' num2str(t) ';']);
else
    evalin('base',['SolverTimes(end+1) = ' num2str(t) ';']);

end

