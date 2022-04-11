global gridN kneeCurrPos kneeCurrVel kneePeakPos
gridN = 2;
kneeCurrPos = 30;
kneeCurrVel = 100;
kneePeakPos = 60;
max_time = 0.3;

init_sim_time = (kneePeakPos - kneeCurrPos)/(kneeCurrVel*gridN);

tic
height_error_min = @(x) (kneePeakPos-x(gridN+1))^2;
% The initial parameter guess; gridN positions, gridN velocities,
% gridN accelerations
x0 = [max_time; linspace(kneeCurrPos,kneePeakPos,gridN)'; linspace(kneeCurrVel,kneeCurrVel,gridN)'; zeros(gridN, 1)];
% No linear inequality or equality constraints
A = [];
b = [];
Aeq = [];
Beq = [];
% Lower bound the position and bound the accelerations between -inf and inf
lb = [0; ones(gridN * 2, 1) * -Inf; ones(gridN, 1) * -1000];
ub = [max_time; ones(gridN * 2, 1) * Inf; ones(gridN, 1) * 1000;];
% Options for fmincon
options = optimoptions(@fmincon, 'TolFun', 1e-3, 'MaxIter', 1e2, ...
                        'ConstraintTolerance', 1e-1,...
                       'MaxFunEvals', 1e4, 'Display', 'iter', ...
                       'Algorithm', 'sqp');
% 'DiffMinChange', 0.001, 
% Solve for the best simulation time + control input
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(height_error_min, x0, A, b, Aeq, Beq, lb, ub, ...
              @double_integrator_constraints,options);
output.message

% Discretize the times
% sim_time = optimal(1);
% delta_time = sim_time / gridN;
% times = 0 : delta_time : sim_time - delta_time;
% Get the state + accelerations (control inputs) out of the vector
sim_time = x(1);
delta_time = sim_time / gridN;
times = 0 : delta_time : sim_time - delta_time;

[positions, vels, accs] = get_traj(x);

% Make the plots
figure(1);
plot(times,accs);
title('Control Input (Acceleration) vs Time');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
figure(2);
plot(times, vels);
title('Velocity vs Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
figure(3);
plot(times, positions);
title('Position vs Time');
xlabel('Time (s)');
ylabel('Position (m)');

disp(sprintf('Finished in %f seconds', toc));

function [ c, ceq ] = double_integrator_constraints( x )
    global gridN kneeCurrPos kneeCurrVel kneePeakPos
    % No nonlinear inequality constraint needed
    c = [];

    % Calculate the timestep
    sim_time = x(1);
    delta_time = sim_time / gridN;

    % Get the states / inputs out of the vector
    [positions, vels, accs] = get_traj(x);
    
    % Constrain initial position and velocity to be zero
    ceq = [positions(1) - kneeCurrPos; vels(1) - kneeCurrVel];
    for i = 1 : length(positions) - 1
        % The state at the beginning of the time interval
        x_i = [positions(i); vels(i)];
        % What the state should be at the start of the next time interval
        x_n = [positions(i+1); vels(i+1)];
        % The time derivative of the state at the beginning of the time
        % interval
        xdot_i = [vels(i); accs(i)];
        % The time derivative of the state at the end of the time interval
        xdot_n = [vels(i+1); accs(i+1)];
        
        % The end state of the time interval calculated using quadrature
        xend = x_i + delta_time * (xdot_i + xdot_n) / 2;
        % Constrain the end state of the current time interval to be
        % equal to the starting state of the next time interval
        ceq = [ceq ; x_n - xend];
    end
    ceq = [ceq ; positions(end) - kneePeakPos; 0];
end

function [positions, vels, accs] = get_traj(x)
global gridN
positions = x(2             : 1 + gridN);
vels      = x(2 + gridN     : 1 + gridN * 2);
accs      = x(2 + gridN * 2 : end);
end