function [x_Hip, y_Hip, x_Knee, y_Knee, x_Foot, y_Foot] = Plot_Position(X, Hip, Knee)

Knee = 180 - Knee;

a = 1; % Second section of leg
b = 1; % First section of leg
c = sqrt(a .^ 2 + b .^ 2 - (2 .* a .* b .* cosd(Knee)));
alpha = asind(a ./ c .* sind(Knee));

x_Hip = zeros(1, length(Hip));
y_Hip = zeros(1, length(Hip));
x_Knee = cosd(180 - Hip) .* b;
y_Knee = sind(180 - Hip) .* b;
x_Foot = cosd(180 - Hip - alpha) .* c;
y_Foot = sind(180 - Hip - alpha) .* c;

% leg = animatedline('MaximumNumPoints', 3);
% leg.LineWidth = 3;
% hipLine = animatedline('MaximumNumPoints', 1);
% hipLine.Marker = 'o';
% kneeLine = animatedline('MaximumNumPoints', 1);
% kneeLine.Marker = 'o';

xhip = zeros(length(X), 1);
yhip = zeros(length(X), 1);

xknee = cosd(180-Hip);
yknee = sind(180-Hip);

q1 = 180 - Hip;
q2 = 180 - Knee;
psi = q1- q2;

xfoot = xknee + cosd(psi);
yfoot = yknee + sind(psi);

% title('Infant Joint Positions over Time');
% 
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% 
% xlim([0, 2]);
% ylim([0, 2]);
% axis square
% set(gcf, 'Position',  [200, 200, 800, 400]) % This changes the location
% of the figure

% Plot animation

% for j = 1:length(xhip)
%     addpoints(leg, xhip(j), yhip(j));
%     addpoints(leg, xknee(j), yknee(j));
%     addpoints(leg, xfoot(j), yfoot(j));
%     
%     addpoints(hipLine, xknee(j), yknee(j));
%     addpoints(kneeLine, xfoot(j), yfoot(j));
%     
%     drawnow;
%     pause(0.05);
%     
% end


end