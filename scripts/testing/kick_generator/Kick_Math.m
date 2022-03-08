%% Generates data for one randomly generated kick




vals = [174.0039    0.0251    0.1248  124.7606    0.0463    2.0604   16.2617    0.0881    2.6480    2.6931    0.2450   -3.7469    1.3434    0.2756   -2.3161    0.7008    0.3128   -1.0626    0.1503    0.4871   -0.8925    5.1932    0.1600   -1.2477];
prep = rand(1, 24) ./ 20 + 1;
vals = vals .* prep;

a1 = vals(1);
b1 = vals(2);
c1 = vals(3);
a2 = vals(4);
b2 = vals(5);
c2 = vals(6);
a3 = vals(7);
b3 = vals(8);
c3 = vals(9);
a4 = vals(10);
b4 = vals(11);
c4 = vals(12);
a5 = vals(13);
b5 = vals(14);
c5 = vals(15);
a6 = vals(16);
b6 = vals(17);
c6 = vals(18);
a7 = vals(19);
b7 = vals(20);
c7 = vals(21);
a8 = vals(22);
b8 = vals(23);
c8 = vals(24);


x = T_Hip;
x = linspace(0,120);
y = a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3) + a4*sin(b4*x+c4) + a5*sin(b5*x+c5) + a6*sin(b6*x+c6) + a7*sin(b7*x+c7) + a8*sin(b8*x+c8);

plot(x,y);
