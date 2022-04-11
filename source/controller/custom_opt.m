% x_0 = 30;
% v_0 = 100;
% x_f = 60;
% v_f = 0;
% T = 0.3;

syms x(t) x_0 v_0 x_f v_f T

alpha = v_f+2*v_0+2/T*(x_0-x_f);
x(t) = (1-(t/T)^2)*(x_0+v_0*t)+(t/T)^2*(x_f+alpha*(t-T))
x(t) = simplify(expand(x(t)))

%%
t_vec = linspace(0,T,100);
x_vec = subs(x(t),t,t_vec);

figure(1)
cla, hold on
plot(t_vec,x_vec)

function plot_tangent_vec(P,P_du)
quiver(P(1),P(2),P_du(1),P_du(2));
end

function new_hermite = hermite_2d(P_0,P_1,P_0_du,P_1_du)
syms new_hermite(u)
new_hermite = [(2*u.^3-3*u.^2+1)*P_0+...
               (-2*u.^3+3*u.^2)*P_1+...
               (u.^3-2*u.^2+u)*P_0_du+...
               (u.^3-u.^2)*P_1_du];
end

function new_bezier = bezier_curve(P)
syms u f(u)
f(u) = 0;
n = size(P,2)-1;
for i = 0:n
    B_i = factorial(n)/(factorial(i)*factorial(n-i))*u^i*(1-u)^(n-i);
    f(u) = f(u) + B_i*P(:,i+1);
end
new_bezier = f(u);
end

function C = format_hermite(H)
H_du = diff(H);
H_ddu = (1/2)*diff(diff(H));
H_dddu = (1/6)*diff(diff(diff(H)));
C_0 = H(0);
C_1 = H_du(0);
C_2 = H_ddu(0);
C_3 = H_dddu(0);
C = [C_3 C_2 C_1 C_0];
end

function plot_2d_curve(curve,num_points)
u_vec = linspace(0,1,num_points);
curve_vec = curve(u_vec);
x = double(curve_vec{1,:});
y = double(curve_vec{2,:});
plot(x,y);
end