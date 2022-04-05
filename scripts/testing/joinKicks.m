function joinKicks
% open the directory where the kicks are
st = dir;
i=3;
newHip = [];
newKnee = [];
calf = 10;
thigh = 10;

while i <= length(st)
    if ~contains(st(i).name, '.csv')
        i = i+1;
    end
    fn = st(i).name;
    
    raw = readmatrix(fn);
    time = [1:length(raw)]';
    hip = raw(:, 3);
    knee = raw(:, 2);
    
%     linspace?? need to figure out how to put leg down and up
% keep same lengths
% circle??
end
end