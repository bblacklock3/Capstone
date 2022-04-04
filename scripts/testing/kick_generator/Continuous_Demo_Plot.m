function [X, Hip, Knee, End_Slope_Hip, End_Slope_Knee, T] = Continuous_Demo_Plot(X, Hip, Knee, Length_Bounds, smooth_Bounds, Diff_Bounds, Angle_End_Bounds, P, isflat, End_Slope_Hip, End_Slope_Knee, h1, h2, leg, visual, T)

% Initializing variables
X_Start = X(end);
Length = rand(1).* (Length_Bounds(2) - Length_Bounds(1)) + Length_Bounds(1);
smooth = rand(1).* (smooth_Bounds(2) - smooth_Bounds(1)) + smooth_Bounds(1);

% Hip joint
Joint = 1;
Angle_Start = Hip(end);
Angle_Diff = rand(1).* (Diff_Bounds(2) - Diff_Bounds(1)) + Diff_Bounds(1);
Angle_End = rand(1).* (Angle_End_Bounds(2) - Angle_End_Bounds(1)) + Angle_End_Bounds(1);
[~, Hip, End_Slope_Hip] = Continuous_Demo_Math(P, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, End_Slope_Hip);

% Knee joint
Joint = 2;
Angle_Start = Knee(end);
Angle_Diff = rand(1).* (Diff_Bounds(4) - Diff_Bounds(3)) + Diff_Bounds(3);
Angle_End = rand(1).* (Angle_End_Bounds(4) - Angle_End_Bounds(3)) + Angle_End_Bounds(3);
[X, Knee, End_Slope_Knee] = Continuous_Demo_Math(P, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, End_Slope_Knee);

% Position plot
[xhip, yhip, xknee, yknee, xfoot, yfoot] = Plot_Position(X, Hip, Knee);

% Plotting flat area
if visual == true
    for i = 1:length(X)
        
        % Angle plot
        addpoints(h1, X(i), Hip(i));
        addpoints(h2, X(i), Knee(i));
        drawnow;
        subplot(1,2,1);
        xlim([X(i) - 3,X(i) + 1]);
        
        % Position plot
        addpoints(leg, xhip(i), yhip(i));
        addpoints(leg, xknee(i), yknee(i));
        addpoints(leg, xfoot(i), yfoot(i));
        
    end
end

% Storing values for final plots
Diff_Hip = diff(Hip) ./ diff(X);
Hip_Diff_Total = [Diff_Hip, Diff_Hip(end)]';
Diff_Knee = diff(Knee) ./ diff(X);
Knee_Diff_Total = [Diff_Knee, Diff_Knee(end)]';
X_Total = X';
Hip_Total = Hip';
Knee_Total = Knee';
Tnew = table(X_Total, Hip_Total, Knee_Total, Hip_Diff_Total, Knee_Diff_Total);
T = [T; Tnew];

end