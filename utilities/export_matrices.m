function [mass, centripetal, coriolis, gravity] = export_matrices(EOM,q)
% Converts the equations of motion into the mass, centripetal, coriolis,
% and gravity matrices.
% Inputs:
% EOM - Column vector of equations of motion
% q   - Variable matrix where each row contains a generalized coordinate
% and each column is the ith derivative with respect to time

EOM_terms = children(EOM);
for i = 1:length(EOM_terms)
    eqns(i) = [EOM_terms{i}{1}];
end
coriolis_vars = nchoosek(q(:,2),2);
coriolis_vars = coriolis_vars(:,1).*coriolis_vars(:,2);
% mass = [];
% centripetal = [];
% coriolis = [];
% gravity = [];
for i = 1:length(eqns)
    for j = 1:size(q,1)
        mass(i,j) = get_coeff(eqns(i),q(j,3));
        centripetal(i,j) = get_coeff(eqns(i),q(j,2)*q(j,2));
    end
    for j = 1:length(coriolis_vars)
        coriolis(i,j) = get_coeff(eqns(i),coriolis_vars(j));
    end
    gravity(i,1) = subs(eqns(i),q(:,2:3),zeros([size(q,1) 2]));
end
mass = simplify(mass);
centripetal = simplify(centripetal);
coriolis = simplify(coriolis);
gravity = simplify(gravity);