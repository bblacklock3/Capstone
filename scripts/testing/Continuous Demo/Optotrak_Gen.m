%% Generates fit values for optotrak kick data sets
% Open Dr Chen's angle data open as the sole open folder

% Resetting MATLAB
clear, clc, close all;

% Changeable variables
filename = 'Fits.mat'; % Filename of output mat file

% Initializing variables
st = dir;
i=3;
type = 'fourier8';
P = {};

% Parsing through each file
while i <= length(st)
    % Checking validity of each file
    if ~contains(st(i).name, '.csv')
        i = i+1;
    else
        
        % Preparing data
        fn = st(i).name;
        Raw = readmatrix(fn);
        Time = [1:length(Raw)]';
        Hip = Raw(:, 3);
        Knee = Raw(:, 2);
        
        % Plotting data
        plot(Time, Hip, Time, Knee);
        
        % Labeling plot
        xlabel('Frames (100/sec)');
        ylabel('Joint Angle (Degrees)');
        title('Infant Joint Angles vs Time');
        lgd = legend('Hip Data', 'Knee Data');
        xlim([min(Time), max(Time)]);
        
        % Prompt for curation of kicks
        prompt = 'Type y for accept, n for skip \n';
        value = input(prompt, 's');
        
        % If kick is selected
        if value == 'y'
            
            % Setting upper and lower bounds for kick
            disp('Select lower bounds');
            [x1, ~] = ginput(1);
            disp('Select upper bounds');
            [x2, ~] = ginput(1);
            hold on;
            xline(x1);
            xline(x2);
            lgd = legend('Hip Data', 'Knee Data');
            hold off;
            
            % Normalizing time range
            idx = find(Time >= x1, 1);
            Time = Time - Time(idx);
            idx = find(Time >= x2 - x1, 1);
            Time = Time ./ Time(idx);
            
            %             % Debugging
            %             figure(2);
            %             plot(Time, Hip, Time, Knee);
            %             hold on;
            %             xline(0);
            %             xline(1);
            
            % Slightly extending fit range
            diff = 0.1;
            x1_p = 0 - diff;
            x2_p = 1 + diff;
            mask = Time >= x1_p & Time <= x2_p;
            Time = Time(mask);
            Hip = Hip(mask);
            Knee = Knee(mask);
            
            % Adjusting fit variable bounds
            polys = numargs(fittype(type)) - 1;
            Bounds = 100;
            Lower = ones(1,polys) .* (0-Bounds);
            Upper = ones(1,polys) .* Bounds;
            
            % Fitting
            p_1 = fit(Time, Hip, type, 'Lower', Lower, 'Upper', Upper);
            p_1 = coeffvalues(p_1)';
            p_2 = fit(Time, Knee, type, 'Lower', Lower, 'Upper', Upper);
            p_2 = coeffvalues(p_2)';
            
            %             % Showing fit
            %             x = linspace(0,1);
            %             y1 = p_1(1)*sin(p_1(2)*x+p_1(3)) + p_1(4)*sin(p_1(5)*x+p_1(6)) + p_1(7)*sin(p_1(8)*x+p_1(9)) + p_1(10)*sin(p_1(11)*x+p_1(12)) + p_1(13)*sin(p_1(14)*x+p_1(15)) + p_1(16)*sin(p_1(17)*x+p_1(18)) + p_1(19)*sin(p_1(20)*x+p_1(21)) + p_1(22)*sin(p_1(23)*x+p_1(24));
            %             y2 = p_2(1)*sin(p_2(2)*x+p_2(3)) + p_2(4)*sin(p_2(5)*x+p_2(6)) + p_2(7)*sin(p_2(8)*x+p_2(9)) + p_2(10)*sin(p_2(11)*x+p_2(12)) + p_2(13)*sin(p_2(14)*x+p_2(15)) + p_2(16)*sin(p_2(17)*x+p_2(18)) + p_2(19)*sin(p_2(20)*x+p_2(21)) + p_2(22)*sin(p_2(23)*x+p_2(24));
            %             figure(3);
            %             plot(x,y1,x,y2);
            
            % Creating array of fits
            P = [P, [p_1, p_2]];
            
        end
        
        % Continuing through files
        i = i+1;
    end
end

% Creating mat file
m = matfile(filename, 'Writable', true);
m.P = P;

