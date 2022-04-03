%% LOAD CONFIG VARIABLES
clear
load([pwd '\CONFIG.mat'])
files = dir(BAG_PATH);

%% FIND ALL BAGFILES
n_files = length(files);
bags.name = {};
bags.path = {};
test.name = {};
test.path = {};
b = 0;
for n = 1:n_files
    [~,~,ext] = fileparts([files(n).name]);
    if strcmp(ext,'.bag')
        b = b + 1;
        bags(b).name = {files(n).name};
        tests(b).name = {files(n).name(1:end-4)};
        bags(b).path = {[BAG_PATH '\' files(n).name]};
        tests(b).path = {[DATA_PATH '\' files(n).name(1:end-4)]};
    end
end

%% CHECK FOR TEST DIR
processed_mask = zeros([1,length(tests)]);
for n = 1:length(tests)
    processed_mask(n) = isfolder(tests(n).path);
end

%% CREATE TEST DIR
process_tests = tests(~processed_mask);
process_bags = bags(~processed_mask);
for n = 1:length(process_tests)
    mkdir(process_tests(n).path{1});
end

%% CONVERT MSGS
for b = 1:length(process_bags)
    process_bags(b).bag = rosbag(process_bags(b).path{1});
    for n = 1:length(TOPIC_LIST)
        clear msg_arr
        [msg_check, msg_arr, time_arr] = ConvertMsgs(process_bags(b),TOPIC_LIST{n});
        if msg_check
            save([process_tests(b).path{1} '\' process_tests(b).name{1} '_' TOPIC_NAMES{n}],'msg_arr','time_arr');
            disp(['Bag: ' process_bags(b).name{1} ' | Topic: ' TOPIC_LIST{n} ' saved!']);
        end
    end
    disp(['Bag: ' process_bags(b).name{1} ' Processed!']);
end


%% LOAD RECEIVED TIMES FOR STD_MSGS
% force_time = [];
% for n = 1:height(bag.MessageList)
%     if bag.MessageList.Topic(n) == '/CEELS/joint_state/force0'
%         force_time = [force_time bag.MessageList.Time(n)];
%     end
% end
% hz = 1/mean(diff(force_time));
% boxplot(1./diff(force_time));


%% Conversion Functions

function [msg_check, msg_arr, time] = ConvertMsgs(process_bag, topic)
topicSel = select(process_bag.bag, 'Topic',topic);
msg_check = true;
time = [];
try
    len = process_bag.bag.AvailableTopics.NumMessages(topic);
    try
        msg_arr = cell(1,len);
        for n = 1:len
            msg = readMessages(topicSel,n);
            msg_arr(n) = msg;
        end
        ts = timeseries(topicSel);
        time = ts.Time';
        disp(['Bag: ' process_bag.name{1} ' | Topic: ' topic ' converted!']);
    catch
        error = ['Read messages failed: ' topic];
        warning(['Bag: ' process_bag.name{1} ' | ' error]);
    end
catch
    error = ['Topic not found: ' topic];
    warning(['Bag: ' process_bag.name{1} ' | ' error]);
    msg_check = false;
    msg_arr = [];
end
end
