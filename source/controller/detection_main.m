function detection_main(t, theta_1, theta_2)
load('ContinuousKicks1.mat');
% T = T(1:4000,:);
t = T.X_Total;
t = t - t(1);
theta_1 = T.Hip_Total;
theta_2 = T.Knee_Total;

% Constants
hip_vel_min = 50;
knee_vel_min = 30;
time_min = 0.1;

% Differentiate joint data
theta_dot_1 = diff(theta_1)./diff(t);
theta_dot_2 = diff(theta_2)./diff(t);
t_diff = t(2:end);

% Apply velocity thresholds
hip_mag = abs(theta_dot_1) > hip_vel_min;
hip_dir = theta_dot_1./abs(theta_dot_1);
knee_mag = abs(theta_dot_2) > knee_vel_min;
knee_dir = theta_dot_2./abs(theta_dot_2);

% Calculate flexion and extension masks
sync_flex = (hip_dir == 1).*(knee_dir == 1).*knee_mag.*hip_mag;
sync_ext = (hip_dir == -1).*(knee_dir == -1).*knee_mag.*hip_mag;

% Find consecutive sections of synchronous motion
sync_flex_conseq = find_conseq(sync_flex);
sync_ext_conseq = find_conseq(sync_ext);

% Apply time threshold
sync_flex_time = t_diff(sync_flex_conseq(:,2))-t_diff(sync_flex_conseq(:,1));
flex_time_mask = find(sync_flex_time > time_min);
sync_flex_conseq = sync_flex_conseq(flex_time_mask,:);
sync_flex_idx = convert_pairs(sync_flex_conseq);
sync_ext_time = t_diff(sync_ext_conseq(:,2))-t_diff(sync_ext_conseq(:,1));
ext_time_mask = find(sync_ext_time > time_min);
sync_ext_conseq = sync_ext_conseq(ext_time_mask,:);
sync_ext_idx = convert_pairs(sync_ext_conseq);

%% Single Plot
% figure(1)
% cla, hold on
% plot(t,theta_1,'-r')
% plot(t,theta_2,'-b')
% plot(t(sync_flex_idx),theta_1(sync_flex_idx),'go')
% plot(t(sync_ext_idx),theta_1(sync_ext_idx),'mo')
% plot(t(sync_flex_idx),theta_2(sync_flex_idx),'go')
% plot(t(sync_ext_idx),theta_2(sync_ext_idx),'mo')

%% Sections
figure(1)
cla, hold on
ax(1) = gca;
plot(t,theta_1,'-r')
plot(t,theta_2,'-b')
min_ang = min([theta_1; theta_2]);
max_ang = max([theta_1; theta_2]);
[flex_x, flex_y] = create_patch(t(sync_flex_conseq), min_ang,max_ang);
[ext_x, ext_y] = create_patch(t(sync_ext_conseq),min_ang,max_ang);
patch(flex_x,flex_y,'g','FaceAlpha',0.3,'EdgeColor','none')
patch(ext_x,ext_y,'m','FaceAlpha',0.3,'EdgeColor','none')
axis('tight')

%% Subplots
figure(2)
ax(2) = subplot(2,1,1);
cla, hold on
plot(t,theta_1,'-r')
plot(t(sync_flex_idx),theta_1(sync_flex_idx),'go')
plot(t(sync_ext_idx),theta_1(sync_ext_idx),'mo')

ax(3) = subplot(2,1,2);
cla, hold on
plot(t,theta_2,'-b')
plot(t(sync_flex_idx),theta_2(sync_flex_idx),'go')
plot(t(sync_ext_idx),theta_2(sync_ext_idx),'mo')

linkaxes(ax,'x');

end

function conseq_idx = find_conseq(vec)
vec_diff = diff(vec);
start_idx = find(vec_diff == 1);
end_idx = find(vec_diff == -1);
if start_idx > end_idx
    end_idx = [end_idx; length(vec_diff)];
end
conseq_idx = [start_idx end_idx]+1;
end

function idx = convert_pairs(pairs)
idx = [];
for i = 1:length(pairs)
    idx = [idx pairs(i,1):pairs(i,2)];
end
end

function [x, y] = create_patch(x_pairs, min_y, max_y)
x = [];
y = [];
for i = 1:length(x_pairs)
    x(:,i) = [x_pairs(i,1); x_pairs(i,2); x_pairs(i,2); x_pairs(i,1)];
    y(:,i) = [min_y; min_y; max_y; max_y];
end
end