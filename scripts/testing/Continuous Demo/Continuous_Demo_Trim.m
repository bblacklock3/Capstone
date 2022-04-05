function [Time, Hip, Knee] = Continuous_Demo_Trim(Time, Hip, Knee, Bounds, Bounds2)

% Initializing variables
istrimmed = false;
Bounds = Bounds .* 0.7;
Bounds2 = Bounds2 .* 0.7;

% Main loop
while istrimmed == false
    index = [];
    
    % Creating velocity and acceleration vectors
    Hip_Diff = diff(Hip) ./ diff(Time);
    Hip_Diff = [Hip_Diff; Hip_Diff(end)];
    Hip_Diff_Diff = diff(Hip_Diff) ./ diff(Time);
    Hip_Diff_Diff = [Hip_Diff_Diff; Hip_Diff_Diff(end)];
    Knee_Diff = diff(Knee) ./ diff(Time);
    Knee_Diff = [Knee_Diff; Knee_Diff(end)];
    Knee_Diff_Diff = diff(Knee_Diff) ./ diff(Time);
    Knee_Diff_Diff = [Knee_Diff_Diff; Knee_Diff_Diff(end)];
    
    % Obtaining mask of out of bounds velocities
    for i = 1:length(Time)
        if Hip_Diff(i) < -Bounds || Hip_Diff(i) > Bounds || Hip_Diff_Diff(i) < -Bounds2 || Hip_Diff_Diff(i) > Bounds2
            index = [index; i];
        elseif Knee_Diff(i) < -Bounds || Knee_Diff(i) > Bounds || Knee_Diff_Diff(i) < -Bounds2 || Knee_Diff_Diff(i) > Bounds2
            index = [index; i];
        end
    end
    
    % Deleting out of bounds velocities
    Time(index) = [];
    Hip(index) = [];
    Knee(index) = [];
    
    % Resetting velocity and acceleration vectors
    Hip_Diff = diff(Hip) ./ diff(Time);
    Hip_Diff = [Hip_Diff; Hip_Diff(end)];
    Hip_Diff_Diff = diff(Hip_Diff) ./ diff(Time);
    Hip_Diff_Diff = [Hip_Diff_Diff; Hip_Diff_Diff(end)];
    Knee_Diff = diff(Knee) ./ diff(Time);
    Knee_Diff = [Knee_Diff; Knee_Diff(end)];
    Knee_Diff_Diff = diff(Knee_Diff) ./ diff(Time);
    Knee_Diff_Diff = [Knee_Diff_Diff; Knee_Diff_Diff(end)];
    
    % Breaking loop if no out of bounds velocities are found
    if isempty(index) == 1
        istrimmed = true;
    end
    
end
end