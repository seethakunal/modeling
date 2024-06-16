tp = [0 25]; % time period over which ode calculated
m = 3; % mass
ICs = [-1.75 0]; % y(0) = -1.75; y'(0) = 0
solveroptions = [];

figure(1)
subplot(2,1,1)

SolverTimes=[];
tic
[t45, y45] = ode45(@(t, y) FES(t, y, m), tp, ICs, solveroptions); 
runtime45 = toc;
times45=SolverTimes;
plot(t45, y45(:,1))

subplot(2,1,2)

SolverTimes=[];
tic
[t15s, y15s] = ode15s(@(t, y) FES(t, y, m), tp, ICs, solveroptions); 
runtime15s = toc;
times15s=SolverTimes;
plot(t15s, y15s(:,1))

figure(2)
bar([1 2], [length(times45) 1000*runtime45; length(times15s) 1000*runtime15s])
xlabel('Solver type')
legend('Solution length', 'Solving time (milliseconds)')
set(gca, 'XTickLabel', {'ode45', 'ode15s'})


figure(3)
hold on
plot(times45, '*')
plot(times15s, 'o')
legend('ode45', 'ode15s')
xlabel('iteration number')
ylabel('time value passed into function')
