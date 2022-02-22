function EOM = calc_lagrange(L,q,Q)

% Calculates the equations of motion with full time derivatives
% L - Lagrange equation
% q - Matrix of symbolic variables where each row represents a generalized
% coordiate
% Q - Virtual work terms
% Example:
%
% % Generalized coordinates
% syms theta theta_dot theta_ddot phi phi_dot phi_ddot s s_dot s_ddot
% % Extra variables
% syms m L g M P
% 
% V_G1 = L/2*theta_dot
% V_G2 = [s_dot; (s-L)*phi_dot]
% 
% m_1 = m;
% m_2 = 2*m;
% 
% I_1 = 1/12*m_1*L^2;
% I_2 = 1/12*m_2*(2*L)^2;
% 
% T = 1/2*m_1*sum(V_G1.^2)+1/2*m_2*sum(V_G2.^2)+1/2*I_1*theta_dot^2+1/2*I_2*phi_dot^2
% V = m*g*L/2*sin(theta)+2*m*g*(s-L)*sin(phi)
% 
% Lagrange = T-V;
% q_diff = [theta theta_dot theta_ddot; phi phi_dot phi_ddot; s s_dot s_ddot;];
% Q = [M; 0; P];
% EOM = calc_lagrange(Lagrange,q_diff,Q)

syms t;
num_q = size(q,1);
syms Ft xt [num_q 3]
for i = 1:num_q
    f_t = str2sym(['f',num2str(i),'(t)']);
    df_t = diff(f_t,t);
    ddf_t = diff(df_t,t);
    Ft(i,:) = [f_t df_t ddf_t];
end
syms EOM [num_q 1];
for i = 1:num_q
    dL_dq_dot = diff(L,q(i,2));
    dL_dq_dot = diff(subs(dL_dq_dot,q,Ft),t);
    dL_dq_dot = subs(dL_dq_dot,Ft,q);
    dL_dq = diff(L,q(i,1));
    EOM(i) = dL_dq_dot-dL_dq;
    EOM(i) = simplify(EOM(i));
    EOM(i) = collect(EOM(i),flip(q(i,:)));
    if nargin == 3
        EOM(i) = EOM(i) == Q(i);
    end
end
end