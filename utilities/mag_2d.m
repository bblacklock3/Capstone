function out = mag_2d(vec)
assert(length(vec)==2)
out = simplify(sqrt(vec(1).^2+vec(2).^2));
end