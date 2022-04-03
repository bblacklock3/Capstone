% Takes in normalized data of angles, outputs two sets of poly values, one
% for each half of the data set with specified fit type
function [P] = Curve_Splitting(Time, Angle, type, num)

P = {};
div = 1/num;
polys = numargs(fittype(type)) - 1;
Lower = ones(1,polys) .* -500000;
Upper = ones(1,polys) .* 500000;
for i = 1:num
    mask = (i-1)*div <= Time & Time <= i*div;
    mask = Time > 0;
    Time_1 = Time(mask);
    Time_1 = Time_1 - Time_1(1);
    Angle_1 = Angle(mask);
    p = fit(Time_1, Angle_1, type, 'Lower', Lower, 'Upper', Upper);
    p = coeffvalues(p);
    P = [P, p];
end

end


% mask = Time < 0.5;
% Time_1 = Time(mask);
% Angle_1 = Angle(mask);
% Time_2 = Time(~mask);
% Time_2 = Time_2 - Time_2(1);
% Angle_2 = Angle(~mask);
%
% P_1 = fit(Time_1, Angle_1, type);
% P_2 = fit(Time_2, Angle_2, type);
% P_1 = coeffvalues(P_1);
% P_2 = coeffvalues(P_2);