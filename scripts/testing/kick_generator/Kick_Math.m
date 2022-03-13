%% Generates data for one randomly generated kick
figure(1);
clf

%pause(5);

for i = 1:100
clf
vals = [174.0039    0.0251    0.1248  124.7606    0.0463    2.0604   16.2617    0.0881    2.6480    2.6931    0.2450   -3.7469    1.3434    0.2756   -2.3161    0.7008    0.3128   -1.0626    0.1503    0.4871   -0.8925    5.1932    0.1600   -1.2477];
%vals = [146.103423226214 0.0490873852123405 0.114999652229547 85.3262684551374 0.098174770424681 1.39815103774157 19.7734716941001 0.196349540849362 1.37662082890395 6.07122656121543 0.294524311274043 1.35163597716866 3.8335063586401 0.392699081698724 1.02416287295631 1.73143744572664 0.490873852123405 0.934183376666444 1.61782043913063 0.589048622548086 0.956968742601002 1.32730641738066 0.687223392972767 0.754327406801219];
%vals = [  166.5049    1.9263    0.8605  103.3772    4.0972    3.0365   10.9344   18.1612   -0.5098    4.5522   27.5556   -4.4748    5.5613   27.8193   -1.6661    0.2920   40.1477    5.2672    3.5779   20.0099  -11.0869    0.1540   68.3153   -1.0036];
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


%x = T_Hip;
x = linspace(0,120);
%x = linspace(0,100);
y = a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3) + a4*sin(b4*x+c4) + a5*sin(b5*x+c5) + a6*sin(b6*x+c6) + a7*sin(b7*x+c7) + a8*sin(b8*x+c8);

plot(x,y);
hold on;
pause(0.01);
end