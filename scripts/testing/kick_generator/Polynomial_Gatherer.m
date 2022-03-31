%% Gathers polynomial for selected curve section

% clear
clc, close all

Data = readmatrix('an01#001.abc.csv');
Time = [1:length(Data)]';
Time = Time - Time(1);
Time = Time ./ Time(end);
Hip = Data(:,3);

X = Time';
Y = Hip';
[p] = Polynomial_Gattherer(X, Y)

function [p] = Polynomial_Gattherer(X, Y)
p = {};
X_Total = {X};
Y_Total = {Y};
lines = [];
diff = 0.05;

for i = 1:5
    check = 'N';
    count = 1;
    while strcmp(check, 'Y') == 0 & strcmp(check, '') == 0
        
        % Plotting raw data
        figure(100);
        hold on;
        for ii = 1:i
            plot(X_Total{ii},Y_Total{ii});
        end
        for iii = 1:length(lines)
            xline(lines(iii), '-.');
        end
        
        % Labeling plot
        xlabel('Time (s)');
        ylabel('Joint Angle (Degrees)');
        title('Infant Joint Angles vs Time');
        
        % Input bounds and poly degree
        if i == 1
            disp('Set Lower x Limit');
            [x1, ~] = ginput(1);
            xline(x1, '-.');
        else
            if count == 1
                x1 = x2;
                x1_original = x1;
            else
                x1 = x1_original;
            end
        end
        disp('Set Upper x Limit');
        [x2, ~] = ginput(1);
        xline(x2, '-.');
        pause(2);
        commandwindow;
        n = input('Input Polynomial Degree: ');
        
        % Slightly extending poly bounds
        x_diff = diff .* (x2 - x1);
        x1_p = x1 - x_diff;
        x2_p = x2 + x_diff;
        
        % Trimming data
        mask = X > x1_p & X < x2_p;
        Trimmed_X = X(mask);
        Trimmed_Y = Y(mask);
        
        % Finding best fit
        p_i = polyfit(Trimmed_X, Trimmed_Y, n);
        X_p = linspace(x1,x2);
        Y_p = polyval(p_i, X_p);
        
        % Plotting fit
        plot(X_p, Y_p);
        
        % Checking if data if correct
        figure(100);
        pause(3);
        commandwindow;
        check = input('Y/N? ', 's');
        
        % Closes figure
        close(figure(100));
        count = count + 1;
    end
    % Preparing output
    p = [p; p_i];
    X_Total = [X_Total X_p];
    Y_Total = [Y_Total Y_p];
    lines = [lines x1 x2];
end

% Polynomial plot
figure(101);








end