function out = isHiMag(hipAng, kneeAng)
mat = corrcoeff(hipAng, kneeAng);
pearson = mat(1, 2);
length = 20;
yknee = length/2*sind(180-hipAng);
q1 = 180 - hipAng;
q2 = 180 - kneeAng;
psi = q1- q2;
yfoot = yknee + length/2*sind(psi); 
out = pearson < 0.3 & yfoot > length*sind(30);
end