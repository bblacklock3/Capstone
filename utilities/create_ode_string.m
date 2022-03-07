function out_str = create_ode_string(name,q_vars,tau_vars,q_ddot_1,q_ddot_2)
q_ddot_1_str = string(q_ddot_1);
q_ddot_2_str = string(q_ddot_2);
assert(length(q_vars) == 4)
for i = 1:4
    q_ddot_1_str = replace(q_ddot_1_str,string(q_vars(i)),['y(' num2str(i) ')']);
    q_ddot_2_str = replace(q_ddot_2_str,string(q_vars(i)),['y(' num2str(i) ')']);
end
for i = 1:length(tau_vars)
    q_ddot_1_str = replace(q_ddot_1_str,string(tau_vars(i)),['tau(' num2str(i) ')']);
    q_ddot_2_str = replace(q_ddot_2_str,string(tau_vars(i)),['tau(' num2str(i) ')']);
end
constants = symvar(q_ddot_1+q_ddot_2);
delete_idx = [];
for i = 1:length(constants)
    if any(constants(i)==q_vars) || any(constants(i)==tau_vars)
        delete_idx = [delete_idx i];
    end
end
constants(delete_idx) = [];
constants_str = string(constants);
constants_str = sort(constants_str);
prop_str = "";
for i = 1:length(constants)
    prop_str = [prop_str constants_str(i) ' = prop.' constants_str(i) ';' newline];
end

out_str = ['function dy = ' name '(t,y,tau,prop)' newline...
    prop_str...
    'dy = ...' newline...
    '[y(2);' newline q_ddot_1_str ';' newline...
    'y(4);' newline q_ddot_2_str ';' newline...
    '];' newline...
    'end'];
out_str = join(out_str,'');
end

