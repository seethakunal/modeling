clear all;
clc;

%% part 1
C = 205;
p = 2.7*10^4;
n = 2;
L = 5;

optionalParams(1) = C;
optionalParams(2) = p;
optionalParams(3) = n;
optionalParams(4) = L;

solverOpts = [];

m=0;

xmesh = linspace(0,L,51);
tspan = linspace(0,100,101);

% part a
sol = pdepe(m,@pdefunc,@ICfunc,@BCfunc,xmesh,tspan,solverOpts, optionalParams);

figure(1);
subplot(1,2,1)
surf(xmesh,tspan,sol)
title("PDE numerical solution, 1a")
xlabel("Length (meters)")
ylabel("Time (seconds)")
zlabel("Temp (K)")

% part b
sol = pdepe(m,@pdefunc,@ICfunc1b,@BCfunc,xmesh,tspan,solverOpts, optionalParams);
subplot(1,2,2)
surf(xmesh,tspan,sol)
title("PDE numerical solution, 1b")
xlabel("Length (meters)")
ylabel("Time (seconds)")
zlabel("Temp (K)")

%% part 2

tspan = linspace(0,10,101); % changed period to 10 seconds

% part a
sol = pdepe(m,@pdefunc,@ICfunc2a,@BCfunc2a,xmesh,tspan,solverOpts, optionalParams);

figure(2);
subplot(1,2,1)
surf(xmesh,tspan,sol)
title("PDE numerical solution, 2a")
xlabel("Length (meters)")
ylabel("Time (seconds)")
zlabel("Temp (K)")

% part b
sol = pdepe(m,@pdefunc,@ICfunc2a,@BCfunc2b,xmesh,tspan,solverOpts, optionalParams);

subplot(1,2,2)
surf(xmesh,tspan,sol)
title("PDE numerical solution, 2b")
xlabel("Length (meters)")
ylabel("Time (seconds)")
zlabel("Temp (K)")


%% functions
function [c, f, s] = pdefunc(x, t, T, dTdx, optionalParams)

C = optionalParams(1);
p = optionalParams(2);
n = optionalParams(3);
L = optionalParams(4);

s=0;
m=0;
c=n*p/C;
f=dTdx;


end

function [T0] = ICfunc(x, optionalParams)

% C = optionalParams(1);
% p = optionalParams(2);
% n = optionalParams(3);
% L = optionalParams(4);

T0=0;

end

function [T0] = ICfunc1b(x, optionalParams)

% C = optionalParams(1);
% p = optionalParams(2);
% n = optionalParams(3);
% L = optionalParams(4);

T0=20*abs(x-2.5)-65;

end

function [pL,qL,pR,qR] = BCfunc(xL,TL,xR,TR, t, optionalParams)

C = optionalParams(1);
p = optionalParams(2);
n = optionalParams(3);
L = optionalParams(4);

pL = TL-(50+t);
qR=0;
qL=0;
pR = TR-(95+t);

% pL=xL;
% qR=50+t;
% qL=95+t;
% pR=xR;

end

function [T0] = ICfunc2a(x, optionalParams)

% C = optionalParams(1);
% p = optionalParams(2);
% n = optionalParams(3);
% L = optionalParams(4);

T0=100;

end

function [pL,qL,pR,qR] = BCfunc2a(xL,TL,xR,TR, t, optionalParams)

C = optionalParams(1);
p = optionalParams(2);
n = optionalParams(3);
L = optionalParams(4);

pL = TL;
qR=-.35;
qL=.35;
pR = TR;

% pL=xL;
% qR=50+t;
% qL=95+t;
% pR=xR;

end

function [pL,qL,pR,qR] = BCfunc2b(xL,TL,xR,TR, t, optionalParams)

C = optionalParams(1);
p = optionalParams(2);
n = optionalParams(3);
L = optionalParams(4);

pL = TL;
qR=.15;
qL=.35;
pR = TR;

% pL=xL;
% qR=50+t;
% qL=95+t;
% pR=xR;

end
