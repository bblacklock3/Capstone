function coeffs = get_coeff(expr,var)
syms sub_var
expr = expand(expr);
sub_expr = subs(expr,var,sub_var);
coeffs = diff(sub_expr,sub_var);
coeff_terms = children(coeffs);
coeff_terms = [coeff_terms{:}];
if any(coeff_terms == sub_var)
    coeffs = 0;
    disp(['Invalid variable for coefficient export'])
end
end