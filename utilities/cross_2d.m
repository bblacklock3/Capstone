function out = cross_2d(vec1,vec2)
switch(length(vec1))
    case 1
        vec1 = [0; 0; vec1];
    case 2
        vec1 = [vec1; 0];
end
switch(length(vec2))
    case 1
        vec2 = [0; 0; vec2];
    case 2
        vec2 = [vec2; 0];
end
out = cross(vec1,vec2);
out = out(1:2);
end

function out = mag_2d(vec)
assert(length(vec)==2)
out = simplify(sqrt(vec(1).^2+vec(2).^2));
end