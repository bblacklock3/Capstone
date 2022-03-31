function plotAllKicks
% cd('C:\..\data\optotrak_data\Chen''s Data\Angles')
st = dir;
% st = dir('Angles');
i=3;
ca = [];
while i <= length(st)
    fn = st(i).name;
    
    % plot the kick
    raw = readmatrix(fn);
    time = [1:length(raw)]';
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
    
    % prompt for the kick type
    prompt = "Type a for async, s for sync, n for not kick, q for quit";
    value = input(prompt);
    ca = [ca; {value}];
    if value ~= 'q'
        i = i + 1;
    else
        break
    end
end

writecell(ca, 'kick_classifications_1.xlsx');
end