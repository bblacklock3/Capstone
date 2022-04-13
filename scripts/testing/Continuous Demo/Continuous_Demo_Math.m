function[x, Angle, End_Slope] = Continuous_Demo_Math(P, Joint, X_Start, Length, Angle_Start, Angle_Diff, Angle_End, smooth, isflat, Init_Slope)




% % Setting normalized coords and gathering delta t
% x = linspace(0,1,(2000.*Length));
% xdiff = x(2) - x(1);
% 
% % Gathering dataset
% r = randi([1,length(P)],1,1);
% P_Angle = P{r}(:,Joint);
% Angle = P_Angle(1) + P_Angle(2)*cos(x*P_Angle(end)) + P_Angle(3)*sin(x*P_Angle(end)) + P_Angle(4)*cos(2*x*P_Angle(end)) + P_Angle(5)*sin(2*x*P_Angle(end)) + P_Angle(6)*cos(3*x*P_Angle(end)) + P_Angle(7)*sin(3*x*P_Angle(end)) + P_Angle(8)*cos(4*x*P_Angle(end)) + P_Angle(9)*sin(4*x*P_Angle(end)) + P_Angle(10)*cos(5*x*P_Angle(end)) + P_Angle(11)*sin(5*x*P_Angle(end)) + P_Angle(12)*cos(6*x*P_Angle(end)) + P_Angle(13)*sin(6*x*P_Angle(end)) + P_Angle(14)*cos(7*x*P_Angle(end)) + P_Angle(15)*sin(7*x*P_Angle(end)) + P_Angle(16)*cos(8*x*P_Angle(end)) + P_Angle(17)*sin(8*x*P_Angle(end));
% Angle_Slope = gradient(Angle) ./ gradient(x);
% Velocity = gradient(Angle) ./ gradient(x);
% Acceleration_Original = gradient(Velocity) ./ gradient(x);
% 
% Velocity_Original = Velocity;
% Angle_Original = Angle;
% x_Original = x;
% xx = x;
% 
% % figure(100);
% % plot(x,Acceleration);
% 
% % Smoothing beginning of angles
% warning off;
% Offset = Velocity(1) .* -1 .* 0.5 .* smooth;
% x_smooth = linspace(-smooth,0,200);
% x_p = linspace(-smooth-0.05,-smooth,100);
% Angle_p = ones(1, length(x_p)) .* (Angle(1) + Offset);
% x_p = [x_p, x(3:end/10)];
% Angle_p = [Angle_p, Angle(3:end/10)];
% p = csapi(x_p, Angle_p, x_smooth);
% warning on;
% x_smooth = x_smooth(1:end-1);
% p = p(1:end-1);
% x = [x_smooth x(3:end)];
% Angle = [p Angle(3:end)];
% Velocity = gradient(Angle) ./ gradient(x);
% Acceleration = gradient(Velocity) ./ gradient(x);
% 
% % Debug plot
% figure(3);
% %plot(x_smooth, p);
% plot(x, Angle);
% 
% figure(4);
% hold on;
% plot(x, Acceleration ./ 10);
% 
% 
% 
% 
% % Adjusting fit variable bounds
% type = 'fourier8';
% polys = numargs(fittype(type)) - 1;
% Bounds = 100;
% Lower = ones(1,polys) .* (0-Bounds);
% Upper = ones(1,polys) .* Bounds;
% 
% flat_x = linspace(-0.1,0,10);
% flat_x = flat_x - smooth;
% flat_Angle = ones(1,10) .* Angle(1);
% x_new = [flat_x, x];
% Angle_new = [flat_Angle, Angle];
% 
% 
% 
% 
% p_1 = fit(x_new', Angle_new', type, 'Lower', Lower, 'Upper', Upper);
% p_1 = coeffvalues(p_1)';
% 
% P_Angle = p_1;
% x = x_new;
% 
% Angle_new = P_Angle(1) + P_Angle(2)*cos(x*P_Angle(end)) + P_Angle(3)*sin(x*P_Angle(end)) + P_Angle(4)*cos(2*x*P_Angle(end)) + P_Angle(5)*sin(2*x*P_Angle(end)) + P_Angle(6)*cos(3*x*P_Angle(end)) + P_Angle(7)*sin(3*x*P_Angle(end)) + P_Angle(8)*cos(4*x*P_Angle(end)) + P_Angle(9)*sin(4*x*P_Angle(end)) + P_Angle(10)*cos(5*x*P_Angle(end)) + P_Angle(11)*sin(5*x*P_Angle(end)) + P_Angle(12)*cos(6*x*P_Angle(end)) + P_Angle(13)*sin(6*x*P_Angle(end)) + P_Angle(14)*cos(7*x*P_Angle(end)) + P_Angle(15)*sin(7*x*P_Angle(end)) + P_Angle(16)*cos(8*x*P_Angle(end)) + P_Angle(17)*sin(8*x*P_Angle(end));
% 
% mask = x_new >= -0.1;
% x_new = x_new(mask);
% Angle_new = Angle_new(mask);
% mask = x_new < -0.01 | x_new > 0.01;
% x_new = x_new(mask);
% Angle_new = Angle_new(mask);
% 
% 
% Velocity_new = gradient(Angle_new) ./ gradient(x_new);
% Acceleration_new = gradient(Velocity_new) ./ gradient(x_new);
% 
% figure(5);
% plot(x_new,Angle_new);
% 
% figure(6);
% plot(x_new, Velocity_new);
% 
% figure(7);
% plot(x_new, Acceleration_new);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Smoothing beginning of acceleration
% warning off;
% x_smooth = linspace(-smooth,0,200);
% x_smooth = [-smooth:xdiff:0];
% x_p = linspace(-smooth-0.05,-smooth,100);
% Angle_p = zeros(1, length(x_p));
% x_p = [x_p, xx(3:end/10)];
% Angle_p = [Angle_p, Acceleration_Original(3:end/10)];
% p = csapi(x_p, Angle_p, x_smooth);
% warning on;
% x_smooth = x_smooth(1:end-1);
% p = p(1:end-1);
% Acceleration2 = [p Acceleration_Original(3:end)];
% 
% xx = [x_smooth xx(3:end)];
% 
% figure(4);
% hold on;
% plot(xx, Acceleration2 ./ 10);
% 
% 
% Vel = cumtrapz(Acceleration2) .* xdiff;
% Pos = cumtrapz(Vel) .* xdiff;
% 
% %Factor = Angle(1) ./ abs(min(Pos));
% %Pos = Pos .* Factor;
% 
% %plot(xx, Vel);
% figure(5);
% hold on;
% plot(xx, Pos);
% 
% figure(6);
% plot(xx, Vel);
% 
% 
% 
% u = sym('u');
% P1 = [-smooth,Angle(1) + Offset; x_Original(1),Angle_Original(1);].';
% 
% P11_du = [0.1,0].';
% P12_du = [1.0,Velocity_Original(3)].';
% 
% H1(u) = cubic_hermite_2d(P1(:,1),P1(:,2),P11_du,P12_du);
% 
% H1_ddu = diff(diff(H1,u),u);
% 
% H1(u) = cubic_hermite_2d(P1(:,1),P1(:,2),P11_du,P12_du);
% 
% figure(7);
% cla, hold on
% plot_2d_curve(H1,100)
% grid on
% legend('Hermite 1')
% 
% 
% %figure(8);
% 
% 
% 
% 
% 
% % P1 = [-smooth,0; 3,8;].';
% % P2 = [3,8; 6,4;].';
% % P3 = [6,4; 11,5;].';
% % P4 = [11,5; x_Original(3),Angle_Original(3)].';
% %
% % P11_du = [1, 0].';
% % P42_du = [1, Velocity_Original(3)].';
% %
% %
% % u = sym('u');
% % P1 = [1,6; 3,8;].';
% % P2 = [3,8; 6,4;].';
% % P3 = [6,4; 11,5;].';
% % P4 = [11,5; 14,1].';
% %
% % P11_du = [1,1].';
% % P42_du = [1,-1].';
% % P12_du = sym('P12_du', [2 1]);
% % P21_du = sym('P21_du', [2 1]);
% % P22_du = sym('P22_du', [2 1]);
% % P31_du = sym('P31_du', [2 1]);
% % P32_du = sym('P32_du', [2 1]);
% % P41_du = sym('P41_du', [2 1]);
% %
% % H1(u) = cubic_hermite_2d(P1(:,1),P1(:,2),P11_du,P12_du);
% % H2(u) = cubic_hermite_2d(P2(:,1),P2(:,2),P21_du,P22_du);
% % H3(u) = cubic_hermite_2d(P3(:,1),P3(:,2),P31_du,P32_du);
% % H4(u) = cubic_hermite_2d(P4(:,1),P4(:,2),P41_du,P42_du);
% %
% % H1_ddu = diff(diff(H1,u),u);
% % H2_ddu = diff(diff(H2,u),u);
% % H3_ddu = diff(diff(H3,u),u);
% % H4_ddu = diff(diff(H4,u),u);
% %
% % C1_12 = P12_du == P21_du;
% % C1_23 = P22_du == P31_du;
% % C1_34 = P32_du == P41_du;
% %
% % C2_12 = H1_ddu(1) == H2_ddu(0);
% % C2_23 = H2_ddu(1) == H3_ddu(0);
% % C2_34 = H3_ddu(1) == H4_ddu(0);
% %
% % s = solve([C1_12 C1_23 C1_34 C2_12 C2_23 C2_34]);
% % P11_du;
% % P12_du = double([s.P12_du1 s.P12_du2].');
% % P21_du = double([s.P21_du1 s.P21_du2].');
% % P22_du = double([s.P22_du1 s.P22_du2].');
% % P31_du = double([s.P31_du1 s.P31_du2].');
% % P32_du = double([s.P32_du1 s.P32_du2].');
% % P41_du = double([s.P41_du1 s.P41_du2].');
% % P42_du;
% %
% % H1(u) = cubic_hermite_2d(P1(:,1),P1(:,2),P11_du,P12_du);
% % H2(u) = cubic_hermite_2d(P2(:,1),P2(:,2),P21_du,P22_du);
% % H3(u) = cubic_hermite_2d(P3(:,1),P3(:,2),P31_du,P32_du);
% % H4(u) = cubic_hermite_2d(P4(:,1),P4(:,2),P41_du,P42_du);
% %
% % figure
% % cla, hold on
% % plot_2d_curve(H1,100)
% % plot_2d_curve(H2,100)
% % plot_2d_curve(H3,100)
% % plot_2d_curve(H4,100)
% % grid on
% % legend('Hermite 1','Hermite 2','Hermite 3','Hermite 4')
% 
% 
%     function new_hermite = cubic_hermite_2d(P_0,P_1,P_0_du,P_1_du)
%         w = sym('u');
%         new_hermite = [(2*w.^3-3*w.^2+1)*P_0+...
%             (-2*w.^3+3*w.^2)*P_1+...
%             (w.^3-2*w.^2+w)*P_0_du+...
%             (w.^3-w.^2)*P_1_du];
%     end




% Setting normalized coords and gathering delta t
x = linspace(0,1,(100.*Length));
xdiff = x(2) - x(1);

% Gathering dataset
r = randi([1,length(P)],1,1);
P_Angle = P{r}(:,Joint);
Angle = P_Angle(1) + P_Angle(2)*cos(x*P_Angle(end)) + P_Angle(3)*sin(x*P_Angle(end)) + P_Angle(4)*cos(2*x*P_Angle(end)) + P_Angle(5)*sin(2*x*P_Angle(end)) + P_Angle(6)*cos(3*x*P_Angle(end)) + P_Angle(7)*sin(3*x*P_Angle(end)) + P_Angle(8)*cos(4*x*P_Angle(end)) + P_Angle(9)*sin(4*x*P_Angle(end)) + P_Angle(10)*cos(5*x*P_Angle(end)) + P_Angle(11)*sin(5*x*P_Angle(end)) + P_Angle(12)*cos(6*x*P_Angle(end)) + P_Angle(13)*sin(6*x*P_Angle(end)) + P_Angle(14)*cos(7*x*P_Angle(end)) + P_Angle(15)*sin(7*x*P_Angle(end)) + P_Angle(16)*cos(8*x*P_Angle(end)) + P_Angle(17)*sin(8*x*P_Angle(end));
Angle_Slope = diff(Angle) ./ diff(x);

% Smoothing beginning of angles
x_smooth = [-smooth:xdiff:0];
n = length(x_smooth);
Angle_Slope_smooth = [Init_Slope:(Angle_Slope(1) - Init_Slope)/n:Angle_Slope(1)];
y_smooth = 0;
Y_smooth = [];
for i = 1:length(Angle_Slope_smooth)
    y_smooth = y_smooth + (Angle_Slope_smooth(i) .* (xdiff));
    Y_smooth = [Y_smooth, y_smooth];
end
Y_smooth = Y_smooth(1:end-1);
Y_smooth = Y_smooth + Angle(1) - Y_smooth(end) + (Angle_Slope(1) .* x_smooth(end));
warning off;
x_p = [x_smooth, x(2:end/10)];
Angle_p = [Y_smooth, Angle(2:end/10)];
warning on;
p = pchip(x_p, Angle_p, x_p);
n = abs(n);
x_p = x_p(1:n);
p = p(1:n);
x = [-smooth(2:end), x_p, x];
Angle = [p, Angle];

% Smoothing end of angles
x_smooth = [1+xdiff:xdiff:1+smooth];
n = length(x_smooth);
Angle_Slope_smooth = [Angle_Slope(end):-Angle_Slope(end)/n:0];
y_smooth = Angle(end);
Y_smooth = [];
for i = 1:length(Angle_Slope_smooth)
    y_smooth = y_smooth + (Angle_Slope_smooth(i) .* (xdiff));
    Y_smooth = [Y_smooth, y_smooth];
end
Y_smooth = Y_smooth(1:end-1);
warning off;
x_p = [x(9*end/10:end), x_smooth];
Angle_p = [Angle(9*end/10:end), Y_smooth];
warning on;
p = pchip(x_p,Angle_p,x_p);
x_p = x_p(end-n+1:end);
p = p(end-n+1:end);
x = [x, x_p(1:end-1)];
Angle = [Angle, p(1:end-1)];

% Formatting angles
if isflat == true
    % Conforming flat angles to input parameters
    Angle = Angle - Angle(1);
    Angle_Diff_T = Angle_Diff;
    Angle_Diff_A = max(abs(max(Angle) - Angle(1)), abs(min(Angle) - Angle(1)));
    Offset = Angle(length(x_smooth)+1);
    Angle(length(x_smooth)+1:end) = Angle(length(x_smooth)+1:end) - Offset;
    Angle(length(x_smooth)+1:end) = Angle(length(x_smooth)+1:end) .* (Angle_Diff_T ./ Angle_Diff_A);
    Angle(length(x_smooth)+1:end) = Angle(length(x_smooth)+1:end) + Offset;
    Angle = Angle + Angle_Start;
    x = x - x(1);
    x = x ./ x(end) .* Length;
    x = x + X_Start;
else
    % Conforming kick angles to input parameters
    Angle = Angle - Angle(1);
    Angle_Diff_T = Angle_Start - Angle_Diff;
    Angle_Diff_A = Angle(1) - min(Angle);
    Angle = Angle .* (Angle_Diff_T ./ Angle_Diff_A);
    Angle = Angle - min(Angle);
    Angle_Diff_T = Angle_End - Angle_Diff;
    Angle_Diff_A = Angle(end) - min(Angle);
    [~,idx] = min(Angle);
    Angle = [Angle(1:idx), Angle(idx+1:end) .* (Angle_Diff_T ./ Angle_Diff_A)];
    Angle = Angle - Angle(1) + Angle_Start;
    x = x - x(1);
    x = x ./ x(end) .* Length;
    x = x + X_Start;
end

% Trimming data and prepping next function call
Angle_Slope = diff(Angle) ./ diff(x);
End_Slope = Angle_Slope(end);
x(1) = [];
Angle(1) = [];

% % Debug plot
% figure(1);
% plot(x, Angle);

end





