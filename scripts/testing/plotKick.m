function plotKick(fn)
% fn is the filename of the .csv file
% 'an01#001.abc.csv'
raw = readmatrix(fn);
time = [1:length(raw)]';
% colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]};

% treat thigh and leg lengths as 10 cm

%% Plotting hip and knee angles over time
subplot(1, 2, 1);
hip = raw(:, 3);
knee = raw(:, 2);


plot(time, hip, time, knee);
hipTrack = animatedline('MaximumNumPoints', 1);
hipTrack.Marker = 'o';
% hipTrack.MarkerEdgeColor = colors{1};
kneeTrack = animatedline('MaximumNumPoints', 1);
kneeTrack.Marker = 'o';
% kneeTrack.MarkerEdgeColor = colors{2};


xlabel('Frames (100/sec)');
ylabel('Joint Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Hip Data', 'Knee Data');
xlim([min(time), max(time)]);
axis square


%% Plotting leg kicking
subplot(1, 2, 2);

leg = animatedline('MaximumNumPoints', 3);
leg.LineWidth = 3;
hipLine = animatedline('MaximumNumPoints', 1);
hipLine.Marker = 'o';
kneeLine = animatedline('MaximumNumPoints', 1);
kneeLine.Marker = 'o';

xhip = zeros(length(time), 1);
yhip = zeros(length(time), 1);

xknee = cosd(180-hip);
yknee = sind(180-hip);

q1 = 180 - hip;
q2 = 180 - knee;
psi = q1- q2;

xfoot = xknee + cosd(psi);
yfoot = yknee + sind(psi);

title('Infant Joint Positions over Time');

set(gca,'xtick',[]);
set(gca,'ytick',[]);

xlim([-1, 2]);
ylim([-1, 2]);
axis square
% set(gcf, 'Position',  [200, 200, 800, 400]) % This changes the location
% of the figure

% Plot animation
for i = 1:1000
    for j = 1:length(xhip)
        addpoints(leg, xhip(j), yhip(j));
        addpoints(leg, xknee(j), yknee(j));
        addpoints(leg, xfoot(j), yfoot(j));
        
        addpoints(hipLine, xknee(j), yknee(j));
        addpoints(kneeLine, xfoot(j), yfoot(j));
        
        addpoints(hipTrack, time(j), hip(j));
        addpoints(kneeTrack, time(j), knee(j));
        
        drawnow update;
        pause(0.05);
        
    end
end