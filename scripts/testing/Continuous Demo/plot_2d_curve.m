function plot_2d_curve(curve,num_points,color)
u_vec = linspace(0,1,num_points);
curve_vec = curve(u_vec);
XX = double(curve_vec{1,:});
YY = double(curve_vec{2,:});
if exist('color','var') == 1
    plot(XX,YY,'Color',color)
else
    plot(XX,YY);
end
end