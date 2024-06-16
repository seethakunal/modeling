clear all;
close all;
clc;

%% script
options = optimset('MaxFunEvals', 1000000, 'MaxIter', 1000000, 'TolFun', 1*10^(-12));
dat = load('FESdata.mat');
time = dat.time;
y0 = [dat.armPosition1, dat.armPosition2, dat.armPosition3];
% y0 one at a time, lower and upper bound at end, initial conditions
lowerBounds = [0 0 0];
upperBounds = [10 1 1];
ICs = [0;0];

global modelFES

% arm position 1
armPos1 = [];
resnorm1 = [];
EstimatedParams1 = [];
residuals1 = [];
for i = 1:25
    initCoeffs = [rand(1)*10, rand(1), rand(1)];
    [EstimatedParams, resnorm, residuals] = lsqcurvefit(@lsqFES, initCoeffs, time, y0(:,1), lowerBounds, upperBounds, options, ICs);
    armPos1(:,i) = modelFES(:,1);
    resnorm1(:,i) = resnorm;
    EstimatedParams1(:,i) = EstimatedParams;
    residuals1(:,i) = residuals;
end
[minval, index] = min(resnorm1);
m = EstimatedParams1(1, index);
ActT = EstimatedParams1(2, index);
ActB = EstimatedParams1(3, index);

solveropts = [];
[modeltime, modelfit] = ode45(@FES, time, ICs, solveropts, m, ActT, ActB);

figure(1)
subplot(3,2,1)
scatter(time, y0(:,1))
hold on
plot(modeltime, modelfit(:,1),'k')
title(['m = ' num2str(m), ', ActT = ' num2str(ActT) ', ActB = ' num2str(ActB)]);

subplot(3,2,2)
stem(time,residuals1(:,index))
title(['Average error (resnorm) = ' num2str(resnorm1(index))])

iterations = 1:25;
figure(2)
subplot(3,1,1)
plot(iterations, resnorm1(:,iterations), 'ro')
title(['Total Avg Error (resnorm) value =' num2str(mean(resnorm1)) '\pm' num2str(std(resnorm1))]);

% armPos2
armPos2 = [];
resnorm2 = [];
EstimatedParams2 = [];
residuals2 = [];
for i = 1:25
    initCoeffs = [rand(1)*10, rand(1), rand(1)];
    [EstimatedParams, resnorm, residuals] = lsqcurvefit(@lsqFES, initCoeffs, time, y0(:,2), lowerBounds, upperBounds, options, ICs);
    armPos2(:,i) = modelFES(:,1);
    resnorm2(:,i) = resnorm;
    EstimatedParams2(:,i) = EstimatedParams;
    residuals2(:,i) = residuals;
end
[minval, index] = min(resnorm2);
m = EstimatedParams2(1, index);
ActT = EstimatedParams2(2, index);
ActB = EstimatedParams2(3, index);

modelfit=[];
solveropts = [];
[modeltime, modelfit] = ode45(@FES, time, ICs, solveropts, m, ActT, ActB);

figure(1)
subplot(3,2,3)
scatter(time, y0(:,2))
hold on
plot(modeltime, modelfit(:,1),'k')
title(['m = ' num2str(m), ', ActT = ' num2str(ActT) ', ActB = ' num2str(ActB)]);

subplot(3,2,4)
stem(time,residuals2(:,index))
title(['Average error (resnorm) = ' num2str(resnorm2(index))])

iterations = 1:25;
figure(2)
subplot(3,1,2)
plot(iterations, resnorm2(:,iterations), 'ro')
title(['Total Avg Error (resnorm) value =' num2str(mean(resnorm2)) '\pm' num2str(std(resnorm2))]);

% armPos3

armPos3 = [];
resnorm3 = [];
EstimatedParams3 = [];
residuals3 = [];
for i = 1:25
    initCoeffs = [rand(1)*10, rand(1), rand(1)];
    [EstimatedParams, resnorm, residuals] = lsqcurvefit(@lsqFES, initCoeffs, time, y0(:,3), lowerBounds, upperBounds, options, ICs);
    armPos3(:,i) = modelFES(:,1);
    resnorm3(:,i) = resnorm;
    EstimatedParams3(:,i) = EstimatedParams;
    residuals3(:,i) = residuals;
end
[minval, index] = min(resnorm3);
m = EstimatedParams3(1, index);
ActT = EstimatedParams3(2, index);
ActB = EstimatedParams3(3, index);

modelfit=[];
solveropts = [];
[modeltime, modelfit] = ode45(@FES, time, ICs, solveropts, m, ActT, ActB);

figure(1)
subplot(3,2,5)
scatter(time, y0(:,3))
hold on
plot(modeltime, modelfit(:,1),'k')
title(['m = ' num2str(m), ', ActT = ' num2str(ActT) ', ActB = ' num2str(ActB)]);

subplot(3,2,6)
stem(time,residuals3(:,index))
title(['Average error (resnorm) = ' num2str(resnorm3(index))])

iterations = 1:25;
figure(2)
subplot(3,1,3)
plot(iterations, resnorm3(:,iterations), 'ro')
title(['Total Avg Error (resnorm) value =' num2str(mean(resnorm3)) '\pm' num2str(std(resnorm3))]);

%% functions

%{
notes
[EstimatedParams, resnorm, residuals] = lsqcurvefit(@ModelFxn,InitCoeffs,xdata,ydata,LB,UB,SolverOptions,ModelParams)
lsqoptions =
optimset('MaxFunEvals',10000000,'MaxIter',100000000000,'TolFun',1E-12);
[coeffs, resnorm, residual] =
lsqcurvefit(@Gompertz,initCoeffs,x,CleanData,lowerBound,upperBound,lsqoptions);
function modelOutput = Gompertz(coeffs,xdata)
This function defines the model of your data
Inputs = experimental input
coeffs = vector of parameters to be estimated
modelOutput = coeffs(1)*exp(coeffs(2)exp...);
for i = 1:10
    initCoeffs(i,:) = lowerBound+(upperBound-lowerBound)
    ....
function modelC = lsqFunc(coeffs,TIME,C0,V,Q)
solverOptions = [];
[t,modelC] = ode45(@odeFun,TIME,c),solverOptions,V,Q,alpha,beta);
end
function dCdt = odeFunc(t,C,V,Q,alpha,beta)
dCdt = (-Q/V)*C = alpha*(C/V)/(beta+C);
end
%}

function [modelFES] = lsqFES(estCoeffs, time, ICs)

global modelFES

solverOpts = [];

% estCoeffs(1) = m, estCoeffs(2) = ActT, estCoeffs(3) = ActB
m = estCoeffs(1);
ActT = estCoeffs(2);
ActB = estCoeffs(3);

% call to ode45 with the estimated parameters
[t, modelFES] = ode45(@FES, time, ICs, solverOpts, m, ActT, ActB);

modelFES = modelFES(:,1);

% save("LSQfitOut.mat", 't', 'estParamSol')

end


function [dy] = FES(t, y, m, ActT, ActB)
% map:
% y = y(1)
% dy(1) = y' = y(2)
% dy(2) = y''

% Force external
Fext = 80*heaviside(t)*exp(-t./2);

% B corresponds to data for biceps, T is for triceps, m is mass
% ActT and ActB given

% resting length
Lrest = 12*(ActB-ActT);

% spring constants
kT = ActT*40;
kB = ActB*20;

% damping constants
cT = 1 + kT/5;
cB = 1 + kB/5;

% calculates the differential equations in the output dy
dy(1) = y(2);
dy(2) = -((cB+cT)/m)*y(2) - ((kB+kT)/m)*y(1) + (kT*Lrest + kB*Lrest)/m + Fext/m;

dy = dy'; % converts dy to column vector by taking transpose

end