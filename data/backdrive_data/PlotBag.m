%% Constants
clear
load('CONFIG.mat')

bags = {};
data_dir = [pwd '\DATA\'];
files = dir(data_dir);
j = 1;
for i = 1:length(files)
    bag_name = files(i).name;
    if length(bag_name) >= 19
        bags(j) = {bag_name};
        j = j+1;
    end
end

%%
bag_name = bags{3};
load([pwd '\DATA\' bag_name '\' bag_name '.mat'])
t_3 = torque.time;
data_3 = torque.data/2;

bag_name = bags{9};
load([pwd '\DATA\' bag_name '\' bag_name '.mat'])
t_9 = torque.time;
data_9 = torque.data/2;

bag_name = bags{10};
load([pwd '\DATA\' bag_name '\' bag_name '.mat'])
t_10 = torque.time;
data_10 = torque.data/2;

bag_name = bags{11};
load([pwd '\DATA\' bag_name '\' bag_name '.mat'])
t_11 = torque.time;
data_11 = torque.data/2;

idx_3 = 3700:3900;
idx_9 = 1250:1450;
idx_10 = 450:650;
idx_11 = 450:650;

t_3 = t_3-t_3(idx_3(1));
t_9 = t_9-t_9(idx_9(1));
t_10 = t_10-t_10(idx_10(1));
t_11 = t_11-t_11(idx_11(1));

figure(1)
cla reset, hold on
plot(t_3(idx_3),data_3(idx_3),'LineWidth',1.5)
plot(t_9(idx_9),data_9(idx_9),'LineWidth',1.5)
plot(t_10(idx_10),data_10(idx_10),'LineWidth',1.5)
plot(t_11(idx_11),data_11(idx_11),'LineWidth',1.5)
legend('Motor Off','I=0 No Calib','Anticogging Calib','Anticogging Calib and Bias');
title('Backdrive Torque')
xlabel('Time (s)')
ylabel('Torque (Nm)')
grid on, grid minor, box on
xlim([0 2])
ylim([0 0.1])

%%
set(0,'DefaultFigureWindowStyle','docked')
for i = 12:length(bags)
    bag_name = bags{i};
    load([pwd '\DATA\' bag_name '\' bag_name '.mat'])
    title_str = 'Backdrive Torque';
    figure(i)
    cla reset, hold on
    plot(torque.time-torque.time(1), torque.data/2)
    xlabel('Time (s)')
    ylabel('Torque (Nm)')
    title(bag_name)
    legend('Measured')
    grid on, grid minor, box on
    % U_s = 2*std(m_torque.data(900:1200)/0.240)
end

%%
bag_name = bags{11};
load([pwd '\DATA\' bag_name '\' bag_name '.mat'])
v = VideoWriter('BackdriveCalibBias','MPEG-4');
% v.FileFormat = 'mp4'
% v.Quality = 100;
v.FrameRate = 30;
open(v);

t = torque.time;
data = torque.data/2;
delta_t = 1/30;
interp_index = 1;
for i = 1:length(t)
    if t(i)-t(interp_index(end)) > delta_t
        interp_index = [interp_index i];
    end
end

font_sz = 20;

figure(11)
cla reset, hold off
for i = 1:length(interp_index)-1
    cla reset, hold on
    plot(t(1:interp_index(i)),data(1:interp_index(i)),'LineWidth',5)
    plot(t(interp_index(i)),data(interp_index(i)),'or','MarkerSize',14,'LineWidth',2)
    xline(t(interp_index(i)),'k','LineWidth',2.5)
    ax = gca;
    ax.FontSize = font_sz;
    title('Backdrive Torque: Anticogging Calibration and Positive Bias','FontSize',font_sz)
    xlabel('Time (s)','FontSize',font_sz)
    ylabel('Torque (Nm)','FontSize',font_sz)
    grid on, grid minor, box on
    xlim([0 10])
    ylim([0 0.02])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
v.Duration
close(v);

%%
bag_name = bags{10};
load([pwd '\DATA\' bag_name '\' bag_name '.mat'])
v = VideoWriter('BackdriveCalib','MPEG-4');
% v.FileFormat = 'mp4'
% v.Quality = 100;
v.FrameRate = 30;
open(v);

t = torque.time;
data = torque.data/2;
delta_t = 1/30;
interp_index = 1;
for i = 1:length(t)
    if t(i)-t(interp_index(end)) > delta_t
        interp_index = [interp_index i];
    end
end

font_sz = 20;

figure(10)
cla reset, hold off
for i = 1:length(interp_index)-1
    cla reset, hold on
    plot(t(1:interp_index(i)),data(1:interp_index(i)),'LineWidth',5)
    plot(t(interp_index(i)),data(interp_index(i)),'or','MarkerSize',14,'LineWidth',2)
    xline(t(interp_index(i)),'k','LineWidth',2.5)
    ax = gca;
    ax.FontSize = font_sz;
    title('Backdrive Torque: Anticogging Calibration','FontSize',font_sz)
    xlabel('Time (s)','FontSize',font_sz)
    ylabel('Torque (Nm)','FontSize',font_sz)
    grid on, grid minor, box on
    xlim([0 10])
    ylim([0 0.032])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
v.Duration
close(v);

%%
bag_name = bags{9};
load([pwd '\DATA\' bag_name '\' bag_name '.mat'])
v = VideoWriter('Backdrive','MPEG-4');
% v.FileFormat = 'mp4'
% v.Quality = 100;
v.FrameRate = 30;
open(v);

t = torque.time;
data = torque.data/2;
delta_t = 1/30;
interp_index = 1;
for i = 1:length(t)
    if t(i)-t(interp_index(end)) > delta_t
        interp_index = [interp_index i];
    end
end

font_sz = 20;

figure(9)
cla reset, hold off
for i = 1:length(interp_index)-1
    cla reset, hold on
    plot(t(1:interp_index(i)),data(1:interp_index(i)),'LineWidth',5)
    plot(t(interp_index(i)),data(interp_index(i)),'or','MarkerSize',14,'LineWidth',2)
    xline(t(interp_index(i)),'k','LineWidth',2.5)
    ax = gca;
    ax.FontSize = font_sz;
    title('Backdrive Torque: Only I_q = 0','FontSize',font_sz)
    xlabel('Time (s)','FontSize',font_sz)
    ylabel('Torque (Nm)','FontSize',font_sz)
    grid on, grid minor, box on
    xlim([0 18])
    ylim([0 0.07])
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
v.Duration
close(v);


%%
% figure(2)
% cla reset, hold on
% [FitModel,Good] = fit(c_torque.data', m_torque.data','poly1');
% plot(FitModel, c_torque.data, m_torque.data)


%%
% figure(3)
% cla reset, hold on
% plot3(c_torque.time, c_torque.data, 1:length(i_torque.time))
% plot3(m_torque.time, m_torque.data, 1:length(m_torque.time))
% view(0,90)

%%
% figure(4)
% cla reset, hold on
% c_seg = 1365:2825;
% m_seg = 1401:2861;
% SaveSeg(m_seg,c_seg,bag_name)
% [FitModel,Good] = fit(c_torque.data(c_seg)', m_torque.data(m_seg)','poly1')
% plot(FitModel, c_torque.data(c_seg), m_torque.data(m_seg))

%%
% figure(5)
% cla reset, hold on
% [FitModel,Good] = fit(c_torque.data(1:end-idx_offset)', m_torque.data(idx_offset:end-1)','poly1')
% plot(c_torque.data(1:end-idx_offset), m_torque.data(idx_offset:end-1),'.','MarkerSize',1)
% xint = linspace(min(c_torque.data(1:end-idx_offset)),max(c_torque.data(1:end-idx_offset)),100);
% CI = predint(FitModel,xint,0.954);
% CI_mean = mean(CI(:,2)-CI(:,1))
% p = plot(xint,FitModel(xint),'-r');
% p(1).LineWidth = 1.5;
% plot(xint,CI,'--','Color',[0.8500 0.3250 0.0980],'MarkerSize',10);
% axis square
% xlim([xint(1) xint(end)])
% ylim([FitModel(xint(1)) FitModel(xint(end))])
% xlabel('Commanded Torque (Nm)')
% ylabel('Measured Torque (Nm)')
% title(title_str)
% 
% L = legend('Data','Fit','K=2','location','southeast');
% % set(L,'visible','off')


function SaveSeg(m_seg,c_seg,bag_name)
save_fn = [pwd '\DATA\' bag_name '\' bag_name '.mat'];
save(save_fn,'m_seg','c_seg','-append')
end
