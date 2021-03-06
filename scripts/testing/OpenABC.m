%% OpenABC

% Resetting MATLAB
clear, clc, close;

% Obtaining data
Angular_Data = xlsread('an01#001.abc.csv'); %Ankle, Knee, Joint
Kinematic_Data = xlsread('LF01#001.ABC.csv'); %Foot, Ankle, Knee, Hip, and Trunk
Time = [1:1:length(Angular_Data)]';
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]};

%% Plot 1

% Plotting
figure(1);
subplot(1,2,1);
hold on;

A1 = Angular_Data(:,1);
A2 = Angular_Data(:,2);
A3 = Angular_Data(:,3);

plot(Time, A1, 'Color', colors{1});
plot(Time, A2, 'Color', colors{2});
plot(Time, A3, 'Color', colors{3});

hA1 = animatedline('MaximumNumPoints', 1);
hA1.Marker = 'o';
hA1.MarkerEdgeColor = colors{1};
hA2 = animatedline('MaximumNumPoints', 1);
hA2.Marker = 'o';
hA2.MarkerEdgeColor = colors{2};
hA3 = animatedline('MaximumNumPoints', 1);
hA3.Marker = 'o';
hA3.MarkerEdgeColor = colors{3};

% Labeling plot
xlabel('Frames (100/sec)');
ylabel('Joint Angle (Degrees)');
title('Infant Joint Angles vs Time');
lgd = legend('Ankle Data', 'Knee Data', 'Hip Data');
xlim([min(Time), max(Time)]);

%% Plot 2

% Obtaining data
figure(1);
subplot(1,2,2);

% Setting plot size
set(gcf, 'Position',  [500, 500, 1000, 400]);

h1 = animatedline('MaximumNumPoints', 6);
h1.LineWidth = 3;
h_Knee = animatedline('MaximumNumPoints', 1);
h_Knee.Marker = 'o';
h_Knee.MarkerEdgeColor = colors{2};
h_Ankle = animatedline('MaximumNumPoints', 1);
h_Ankle.Marker = 'o';
h_Ankle.MarkerEdgeColor = colors{1};
h_Hip = animatedline('MaximumNumPoints', 1);
h_Hip.Marker = 'o';
h_Hip.MarkerEdgeColor = colors{3};

x1 = Kinematic_Data(:,1);
y1 = Kinematic_Data(:,2);
z1 = Kinematic_Data(:,3);

x2 = Kinematic_Data(:,4);
y2 = Kinematic_Data(:,5);
z2 = Kinematic_Data(:,6);

x3 = Kinematic_Data(:,7);
y3 = Kinematic_Data(:,8);
z3 = Kinematic_Data(:,9);

x4 = Kinematic_Data(:,10);
y4 = Kinematic_Data(:,11);
z4 = Kinematic_Data(:,12);

x5 = Kinematic_Data(:,13);
y5 = Kinematic_Data(:,14);
z5 = Kinematic_Data(:,15);

x6 = Kinematic_Data(:,16);
y6 = Kinematic_Data(:,17);
z6 = Kinematic_Data(:,18);

% Plot setup
view(3);
view(-100,45);
xlim([-100, 30]);
ylim([0, 300]);
zlim([-2420, -2280]);

% Label plot
%xlabel('Position (mm)');
%ylabel('Position (mm)');
%zlabel('Position (mm)');
title('Infant Joint Positions over Time');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

% Plot animation
for i = 1:1000
    for k = 1:length(x1)
        addpoints(h1,x1(k),y1(k),z1(k));
        addpoints(h1,x2(k),y2(k),z2(k));
        addpoints(h1,x3(k),y3(k),z3(k));
        addpoints(h1,x4(k),y4(k),z4(k));
        addpoints(h1,x5(k),y5(k),z5(k));
        addpoints(h1,x6(k),y6(k),z6(k));
        
        addpoints(h_Knee,x3(k),y3(k),z3(k));
        addpoints(h_Ankle,x2(k),y2(k),z2(k));
        addpoints(h_Hip,x4(k),y4(k),z4(k));
        
        addpoints(hA1,Time(k),A1(k));
        addpoints(hA2,Time(k),A2(k));
        addpoints(hA3,Time(k),A3(k));
        
        
        drawnow update;
        
        
        
        
        
        pause(0.05);
    end
end




