function ExtractData(bag_name)
%% Constants
load('CONFIG.mat')

%% Load Data
tension = [];
torque = [];

%% Current
disp('Loading Tension')
load(GetTopic('tension',bag_name))
tension.time = time_arr;
tension.data = zeros([1 length(msg_arr)]);
for i = 1:length(msg_arr)
tension.data(i) = msg_arr{i}.Data;
end

disp('Loading Torque')
load(GetTopic('torque',bag_name))
torque.time = time_arr;
torque.data = zeros([1 length(msg_arr)]);
for i = 1:length(msg_arr)
torque.data(i) = msg_arr{i}.Data;
end

%% Offset by start time
start_times = [tension.time(1) torque.time(1)];
start_time = inf;
for i = 1:length(start_time)
    t = start_times(i);
    if t < start_time; start_time = t; end
end
tension.time = tension.time - start_time;
torque.time = torque.time - start_time;

%% Trim vectors
max_len = inf;
lengths = [length(tension.time) length(torque.time)];
for i = 1:length(lengths)
    l = lengths(i);
    if l < max_len; max_len = l; end
end
tension.time = tension.time(1:max_len);
torque.time = torque.time(1:max_len);
tension.data = tension.data(1:max_len);
torque.data = torque.data(1:max_len);

save_fn = [pwd '\DATA\' bag_name '\' bag_name '.mat'];
save(save_fn,'tension','torque')
end

function name = GetTopic(topic, bag_name)
name = [pwd '\DATA\' bag_name '\' bag_name '_' topic '.mat'];
end

