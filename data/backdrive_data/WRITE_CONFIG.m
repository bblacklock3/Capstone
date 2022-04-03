clear
DATA_PATH = [pwd '\DATA'];
BAG_PATH = [DATA_PATH '\BAGS'];
TOPIC_LIST = {'/torque'
    };
TOPIC_NAMES = {};
for n = 1:length(TOPIC_LIST)
    topic = TOPIC_LIST{n};
    TOPIC_NAMES = [TOPIC_NAMES {strrep(topic(2:end),'/','-')}];
end
save('CONFIG.mat','DATA_PATH','BAG_PATH','TOPIC_LIST','TOPIC_NAMES');