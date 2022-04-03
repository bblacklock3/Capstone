function ExtractTorque(bag_name)
%% Constants
load('CONFIG.mat')

%% Load Data
torque = [];

disp('Loading Torque')
load(GetTopic('torque',bag_name))
torque.time = time_arr;
torque.data = zeros([1 length(msg_arr)]);
for i = 1:length(msg_arr)
torque.data(i) = msg_arr{i}.Data;
end

save_fn = [pwd '\DATA\' bag_name '\' bag_name '.mat'];
save(save_fn,'torque')
end

function name = GetTopic(topic, bag_name)
name = [pwd '\DATA\' bag_name '\' bag_name '_' topic '.mat'];
end

