%% Automatically generates plots of n random kicks

% Resetting MATLAB
clear, clc, close all;

load('P.mat');


Dy1 = {};
Dy2 = {};
for i = 1:length(P)
    
    %y1 = P_Hip(1)*sin(P_Hip(2)*x+P_Hip(3)) + P_Hip(4)*sin(P_Hip(5)*x+P_Hip(6)) + P_Hip(7)*sin(P_Hip(8)*x+P_Hip(9)) + P_Hip(10)*sin(P_Hip(11)*x+P_Hip(12)) + P_Hip(13)*sin(P_Hip(14)*x+P_Hip(15)) + P_Hip(16)*sin(P_Hip(17)*x+P_Hip(18)) + P_Hip(19)*sin(P_Hip(20)*x+P_Hip(21)) + P_Hip(22)*sin(P_Hip(23)*x+P_Hip(24));
    %y2 = P_Knee(1)*sin(P_Knee(2)*x+P_Knee(3)) + P_Knee(4)*sin(P_Knee(5)*x+P_Knee(6)) + P_Knee(7)*sin(P_Knee(8)*x+P_Knee(9)) + P_Knee(10)*sin(P_Knee(11)*x+P_Knee(12)) + P_Knee(13)*sin(P_Knee(14)*x+P_Knee(15)) + P_Knee(16)*sin(P_Knee(17)*x+P_Knee(18)) + P_Knee(19)*sin(P_Knee(20)*x+P_Knee(21)) + P_Knee(22)*sin(P_Knee(23)*x+P_Knee(24));
    %y1 = P_Hip(1) + P_Hip(2)*cos(x*P_Hip(end)) + P_Hip(3)*sin(x*P_Hip(end)) + P_Hip(4)*cos(2*x*P_Hip(end)) + P_Hip(5)*sin(2*x*P_Hip(end)) + P_Hip(6)*cos(3*x*P_Hip(end)) + P_Hip(7)*sin(3*x*P_Hip(end)) + P_Hip(8)*cos(4*x*P_Hip(end)) + P_Hip(9)*sin(4*x*P_Hip(end)) + P_Hip(10)*cos(5*x*P_Hip(end)) + P_Hip(11)*sin(5*x*P_Hip(end)) + P_Hip(12)*cos(6*x*P_Hip(end)) + P_Hip(13)*sin(6*x*P_Hip(end)) + P_Hip(14)*cos(7*x*P_Hip(end)) + P_Hip(15)*sin(7*x*P_Hip(end)) + P_Hip(16)*cos(8*x*P_Hip(end)) + P_Hip(17)*sin(8*x*P_Hip(end));
    %y2 = P_Knee(1) + P_Knee(2)*cos(x*P_Knee(end)) + P_Knee(3)*sin(x*P_Knee(end)) + P_Knee(4)*cos(2*x*P_Knee(end)) + P_Knee(5)*sin(2*x*P_Knee(end)) + P_Knee(6)*cos(3*x*P_Knee(end)) + P_Knee(7)*sin(3*x*P_Knee(end)) + P_Knee(8)*cos(4*x*P_Knee(end)) + P_Knee(9)*sin(4*x*P_Knee(end)) + P_Knee(10)*cos(5*x*P_Knee(end)) + P_Knee(11)*sin(5*x*P_Knee(end)) + P_Knee(12)*cos(6*x*P_Knee(end)) + P_Knee(13)*sin(6*x*P_Knee(end)) + P_Knee(14)*cos(7*x*P_Knee(end)) + P_Knee(15)*sin(7*x*P_Knee(end)) + P_Knee(16)*cos(8*x*P_Knee(end)) + P_Knee(17)*sin(8*x*P_Knee(end));
    
    %y1_ii = P_Hip(1) + P_Hip(2)*cos(x*P_Hip(end)) + P_Hip(3)*sin(x*P_Hip(end)) + P_Hip(4)*cos(2*x*P_Hip(end)) + P_Hip(5)*sin(2*x*P_Hip(end)) + P_Hip(6)*cos(3*x*P_Hip(end)) + P_Hip(7)*sin(3*x*P_Hip(end)) + P_Hip(8)*cos(4*x*P_Hip(end)) + P_Hip(9)*sin(4*x*P_Hip(end)) + P_Hip(10)*cos(5*x*P_Hip(end)) + P_Hip(11)*sin(5*x*P_Hip(end));
    %y2_ii = P_Knee(1) + P_Knee(2)*cos(x*P_Knee(end)) + P_Knee(3)*sin(x*P_Knee(end)) + P_Knee(4)*cos(2*x*P_Knee(end)) + P_Knee(5)*sin(2*x*P_Knee(end)) + P_Knee(6)*cos(3*x*P_Knee(end)) + P_Knee(7)*sin(3*x*P_Knee(end)) + P_Knee(8)*cos(4*x*P_Knee(end)) + P_Knee(9)*sin(4*x*P_Knee(end)) + P_Knee(10)*cos(5*x*P_Knee(end)) + P_Knee(11)*sin(5*x*P_Knee(end));
    
    %y1_ii = P_Hip(1) + P_Hip(2)*cos(x*P_Hip(end)) + P_Hip(3)*sin(x*P_Hip(end)) + P_Hip(4)*cos(2*x*P_Hip(end)) + P_Hip(5)*sin(2*x*P_Hip(end)) + P_Hip(6)*cos(3*x*P_Hip(end)) + P_Hip(7)*sin(3*x*P_Hip(end));
    %y2_ii = P_Knee(1) + P_Knee(2)*cos(x*P_Knee(end)) + P_Knee(3)*sin(x*P_Knee(end)) + P_Knee(4)*cos(2*x*P_Knee(end)) + P_Knee(5)*sin(2*x*P_Knee(end)) + P_Knee(6)*cos(3*x*P_Knee(end)) + P_Knee(7)*sin(3*x*P_Knee(end));
    
    
    %         if ii == 1
    %             y1_offset = y1_ii(1);
    %             y2_offset = y2_ii(1);
    %         end
    %
    %         y1_ii = y1_ii + y1_diff;
    %         y2_ii = y2_ii + y2_diff;
    %         y1_diff = y1_ii(1) - y1_offset;
    %         y2_diff = y2_ii(1) - y2_offset;
    %
    %         if ii ~= 1
    %             y1_offset = y1_ii(end);
    %             y2_offset = y2_ii(end);
    %         end
    
    
    P_Hip_Original = P{i,1};
    P_Knee_Original = P{i,2};
    y1 = [];
    y2 = [];
    X = [];
    Dy1_ii = {};
    Dy2_ii = {};
    for ii = 1:length(P{i,1})
        P_Hip = P_Hip_Original{ii};
        P_Knee = P_Knee_Original{ii};
        
        x = linspace(0,1/length(P{i,1}));
        y1_ii = P_Hip(1)*sin(P_Hip(2)*x+P_Hip(3)) + P_Hip(4)*sin(P_Hip(5)*x+P_Hip(6)) + P_Hip(7)*sin(P_Hip(8)*x+P_Hip(9)) + P_Hip(10)*sin(P_Hip(11)*x+P_Hip(12)) + P_Hip(13)*sin(P_Hip(14)*x+P_Hip(15)) + P_Hip(16)*sin(P_Hip(17)*x+P_Hip(18)) + P_Hip(19)*sin(P_Hip(20)*x+P_Hip(21)) + P_Hip(22)*sin(P_Hip(23)*x+P_Hip(24));
        y2_ii = P_Knee(1)*sin(P_Knee(2)*x+P_Knee(3)) + P_Knee(4)*sin(P_Knee(5)*x+P_Knee(6)) + P_Knee(7)*sin(P_Knee(8)*x+P_Knee(9)) + P_Knee(10)*sin(P_Knee(11)*x+P_Knee(12)) + P_Knee(13)*sin(P_Knee(14)*x+P_Knee(15)) + P_Knee(16)*sin(P_Knee(17)*x+P_Knee(18)) + P_Knee(19)*sin(P_Knee(20)*x+P_Knee(21)) + P_Knee(22)*sin(P_Knee(23)*x+P_Knee(24));
        dy1_ii = gradient(y1_ii) ./ gradient(x);
        dy2_ii = gradient(y2_ii) ./ gradient(x);
        dy1_ii = [dy1_ii(1), dy1_ii(end)];
        dy2_ii = [dy2_ii(1), dy2_ii(end)];
        
        if ii == 1
            y1_ii_end = y1_ii(1);
            y2_ii_end = y2_ii(1);
        end
        y1_offset = y1_ii(1) - y1_ii_end;
        y2_offset = y2_ii(1) - y2_ii_end;
        y1_ii = y1_ii - y1_offset;
        y2_ii = y2_ii - y2_offset;
        y1_ii_end = y1_ii(end);
        y2_ii_end = y2_ii(end);
        
        y1 = [y1, y1_ii];
        y2 = [y2, y2_ii];
        Dy1_ii = [Dy1_ii, dy1_ii];
        Dy2_ii = [Dy2_ii, dy2_ii];
        x = x + x(end)*(ii-1);
        X = [X, x];
        
    end
    Dy1 = [Dy1; Dy1_ii];
    Dy2 = [Dy2; Dy2_ii];
    
    figure(i);
    hold on;
    plot(X, y1);
    plot(X, y2);
    
    
    for iii = 1:10
        canpass = 0;
        while canpass == 0
            y1 = [];
            y2 = [];
            X = [];
            Dy1_ii = {};
            Dy2_ii = {};
            Per = 10;
            for ii = 1:length(P{i,1})
                pass = 0;
                while pass == 0
                    P_Hip = P_Hip_Original{ii};
                    random = (rand(1, length(P_Hip)) .* Per ./ 100) + 1;
                    P_Hip = P_Hip.*random;
                    P_Knee = P_Knee_Original{ii};
                    random = (rand(1, length(P_Knee)) .* Per ./ 100) + 1;
                    P_Knee = P_Knee.*random;
                    
                    x = linspace(0,1/length(P{i,1}));
                    y1_ii = P_Hip(1)*sin(P_Hip(2)*x+P_Hip(3)) + P_Hip(4)*sin(P_Hip(5)*x+P_Hip(6)) + P_Hip(7)*sin(P_Hip(8)*x+P_Hip(9)) + P_Hip(10)*sin(P_Hip(11)*x+P_Hip(12)) + P_Hip(13)*sin(P_Hip(14)*x+P_Hip(15)) + P_Hip(16)*sin(P_Hip(17)*x+P_Hip(18)) + P_Hip(19)*sin(P_Hip(20)*x+P_Hip(21)) + P_Hip(22)*sin(P_Hip(23)*x+P_Hip(24));
                    y2_ii = P_Knee(1)*sin(P_Knee(2)*x+P_Knee(3)) + P_Knee(4)*sin(P_Knee(5)*x+P_Knee(6)) + P_Knee(7)*sin(P_Knee(8)*x+P_Knee(9)) + P_Knee(10)*sin(P_Knee(11)*x+P_Knee(12)) + P_Knee(13)*sin(P_Knee(14)*x+P_Knee(15)) + P_Knee(16)*sin(P_Knee(17)*x+P_Knee(18)) + P_Knee(19)*sin(P_Knee(20)*x+P_Knee(21)) + P_Knee(22)*sin(P_Knee(23)*x+P_Knee(24));
                    dy1_ii = gradient(y1_ii) ./ gradient(x);
                    dy2_ii = gradient(y2_ii) ./ gradient(x);
                    dy1_ii = [dy1_ii(1), dy1_ii(end)];
                    dy2_ii = [dy2_ii(1), dy2_ii(end)];
                    if ii == 1
                        y1_ii_end = y1_ii(1);
                        y2_ii_end = y2_ii(1);
                    end
                    y1_offset = y1_ii(1) - y1_ii_end;
                    y2_offset = y2_ii(1) - y2_ii_end;
                    y1_ii = y1_ii - y1_offset;
                    y2_ii = y2_ii - y2_offset;
                    
                    dy1 = Dy1{i,ii};
                    dy2 = Dy2{i,ii};
                    dy1_first = dy1(1);
                    dy2_first = dy2(1);
                    dy1_last = dy1(2);
                    dy2_last = dy2(2);
                    
                    %                 if abs(dy1_first) >= (0.9 * abs(dy1_ii(1))) && abs(dy1_first) <= (1.1 * abs(dy1_ii(1)))
                    %                     if abs(dy1_last) >= (0.9 * abs(dy1_ii(2))) && abs(dy1_last) <= (1.1 * abs(dy1_ii(2)))
                    %                         pass = 1;
                    %                         if ii ~= 1
                    %                             % abs(dy1_ii(1))
                    %                             %if abs(Dy1_ii{ii-1}(2)) >= (0.5 * abs(dy1_ii(1))) && abs(Dy1_ii{ii-1}(2)) <= (2 * abs(dy1_ii(1)))
                    %                             %    pass = 1;
                    %                             %else
                    %                             %     pass = 0;
                    %                             %end
                    %                         end
                    %                     else
                    %                         pass = 0;
                    %                     end
                    %                 else
                    %                     pass = 0;
                    %                 end
                    
                    if max(y1_ii) > 170 || max(y2_ii) > 170 || min(y1_ii) < 60 || min(y2_ii) < 60
                        pass = 0;
                        Per = Per .* 0.95;
                    else
                        pass = 1;
                        Per = 5;
                    end
                    
                    
                    
                    %                 dy1_ii
                    %                 dy2_ii
                    %
                    %                 if dy1_ii(1) > -0.1 && dy1_ii(1) < 0.1
                    %                     pass = 0;
                    %                 elseif dy1_ii(2) > -0.1 && dy1_ii(2) < 0.1
                    %                     pass = 0;
                    %                 elseif dy2_ii(1) > -0.1 && dy2_ii(1) < 0.1
                    %                     pass = 0;
                    %                 elseif dy2_ii(2) > -0.1 && dy2_ii(2) < 0.1
                    %                     pass = 0;
                    %                 end
                    
                    if pass == 1
                        y1_ii_end = y1_ii(end);
                        y2_ii_end = y2_ii(end);
                        y1 = [y1, y1_ii];
                        y2 = [y2, y2_ii];
                        Dy1_ii = [Dy1_ii, dy1_ii];
                        Dy2_ii = [Dy2_ii, dy2_ii];
                        x = x + x(end)*(ii-1);
                        X = [X, x];
                    end
                    
                end
            end
            
            dT = x(2) - x(1);
            Fs = 1/dT;
            
            X = [linspace(-.2,-0.001) X linspace(1.001,1.2)];
            y1 = [ones(1,100).*y1(1), y1, ones(1,100).*y1(end)];
            y1 = lowpass(y1, 0.001, Fs);
            mask = X >= 0 & X <= 1;
            y1 = y1(mask);
            
            y2 = [ones(1,100).*y2(1), y2, ones(1,100).*y2(end)];
            y2 = lowpass(y2, 0.001, Fs);
            y2 = y2(mask);
            X = X(mask);
            
            dy1 = gradient(y1) ./ gradient(X);
            dy2 = gradient(y2) ./ gradient(X);
            
            if dy1(1) > 100 || dy1(1) < -100
            elseif dy1(end) > 100 || dy1(end) < -100
            elseif dy2(1) > 100 || dy2(1) < -100
            elseif dy2(end) > 100 || dy2(end) < -100
            else
                canpass = 1;
            end
            
            
        end
        
        figure(i);
        hold on;
        plot(X, y1);
        plot(X, y2);
        
    end
    
    
    
    %     P_Hip = P{i,1};
    %     P_Knee = P{i,2};
    %     P_Hip_Original_1 = P_Hip;
    %     P_Knee_Original_1 = P_Knee;
    %     x = linspace(0,0.5);
    %     y1 = P_Hip(1) + P_Hip(2)*cos(x*P_Hip(end)) + P_Hip(3)*sin(x*P_Hip(end)) + P_Hip(4)*cos(2*x*P_Hip(end)) + P_Hip(5)*sin(2*x*P_Hip(end)) + P_Hip(6)*cos(3*x*P_Hip(end)) + P_Hip(7)*sin(3*x*P_Hip(end)) + P_Hip(8)*cos(4*x*P_Hip(end)) + P_Hip(9)*sin(4*x*P_Hip(end)) + P_Hip(10)*cos(5*x*P_Hip(end)) + P_Hip(11)*sin(5*x*P_Hip(end));
    %     y2 = P_Knee(1) + P_Knee(2)*cos(x*P_Knee(end)) + P_Knee(3)*sin(x*P_Knee(end)) + P_Knee(4)*cos(2*x*P_Knee(end)) + P_Knee(5)*sin(2*x*P_Knee(end)) + P_Knee(6)*cos(3*x*P_Knee(end)) + P_Knee(7)*sin(3*x*P_Knee(end)) + P_Knee(8)*cos(4*x*P_Knee(end)) + P_Knee(9)*sin(4*x*P_Knee(end)) + P_Knee(10)*cos(5*x*P_Knee(end)) + P_Knee(11)*sin(5*x*P_Knee(end));
    %
    %     P_Hip = P{i,3};
    %     P_Knee = P{i,4};
    %     P_Hip_Original_2 = P_Hip;
    %     P_Knee_Original_2 = P_Knee;
    %     y3 = P_Hip(1) + P_Hip(2)*cos(x*P_Hip(end)) + P_Hip(3)*sin(x*P_Hip(end)) + P_Hip(4)*cos(2*x*P_Hip(end)) + P_Hip(5)*sin(2*x*P_Hip(end)) + P_Hip(6)*cos(3*x*P_Hip(end)) + P_Hip(7)*sin(3*x*P_Hip(end)) + P_Hip(8)*cos(4*x*P_Hip(end)) + P_Hip(9)*sin(4*x*P_Hip(end)) + P_Hip(10)*cos(5*x*P_Hip(end)) + P_Hip(11)*sin(5*x*P_Hip(end));
    %     y4 = P_Knee(1) + P_Knee(2)*cos(x*P_Knee(end)) + P_Knee(3)*sin(x*P_Knee(end)) + P_Knee(4)*cos(2*x*P_Knee(end)) + P_Knee(5)*sin(2*x*P_Knee(end)) + P_Knee(6)*cos(3*x*P_Knee(end)) + P_Knee(7)*sin(3*x*P_Knee(end)) + P_Knee(8)*cos(4*x*P_Knee(end)) + P_Knee(9)*sin(4*x*P_Knee(end)) + P_Knee(10)*cos(5*x*P_Knee(end)) + P_Knee(11)*sin(5*x*P_Knee(end));
    %     x2 = x + 0.5;
    %     diff = y1(end) - y3(1);
    %     y3 = y3 + diff;
    %     diff = y2(end) - y4(1);
    %     y4 = y4 + diff;
    %
    %     figure(i);
    %     hold on;
    %     plot([x, x2],[y1, y3]);
    %     plot([x, x2],[y2, y4]);
    
    %     for ii = 1:100
    %
    %         random = rand(1, length(P_Hip_Original_1)) ./ 200 + 1;
    %         %random = ones(1, length(P_Hip_Original_1));
    %         P_Hip = P_Hip_Original_1.* random;
    %         P_Knee = P_Knee_Original_1 .* random;
    %         x = linspace(0,0.5);
    %         y1 = P_Hip(1) + P_Hip(2)*cos(x*P_Hip(end)) + P_Hip(3)*sin(x*P_Hip(end)) + P_Hip(4)*cos(2*x*P_Hip(end)) + P_Hip(5)*sin(2*x*P_Hip(end)) + P_Hip(6)*cos(3*x*P_Hip(end)) + P_Hip(7)*sin(3*x*P_Hip(end)) + P_Hip(8)*cos(4*x*P_Hip(end)) + P_Hip(9)*sin(4*x*P_Hip(end)) + P_Hip(10)*cos(5*x*P_Hip(end)) + P_Hip(11)*sin(5*x*P_Hip(end));
    %         y2 = P_Knee(1) + P_Knee(2)*cos(x*P_Knee(end)) + P_Knee(3)*sin(x*P_Knee(end)) + P_Knee(4)*cos(2*x*P_Knee(end)) + P_Knee(5)*sin(2*x*P_Knee(end)) + P_Knee(6)*cos(3*x*P_Knee(end)) + P_Knee(7)*sin(3*x*P_Knee(end)) + P_Knee(8)*cos(4*x*P_Knee(end)) + P_Knee(9)*sin(4*x*P_Knee(end)) + P_Knee(10)*cos(5*x*P_Knee(end)) + P_Knee(11)*sin(5*x*P_Knee(end));
    %
    %         %random = rand(1, length(P_Hip_Original_2)) ./ 200000 + 1;
    %         random = ones(1, length(P_Hip_Original_2));
    %         P_Hip = P_Hip_Original_2.* random;
    %         P_Knee = P_Knee_Original_2 .* random;
    %         y3 = P_Hip(1) + P_Hip(2)*cos(x*P_Hip(end)) + P_Hip(3)*sin(x*P_Hip(end)) + P_Hip(4)*cos(2*x*P_Hip(end)) + P_Hip(5)*sin(2*x*P_Hip(end)) + P_Hip(6)*cos(3*x*P_Hip(end)) + P_Hip(7)*sin(3*x*P_Hip(end)) + P_Hip(8)*cos(4*x*P_Hip(end)) + P_Hip(9)*sin(4*x*P_Hip(end)) + P_Hip(10)*cos(5*x*P_Hip(end)) + P_Hip(11)*sin(5*x*P_Hip(end));
    %         y4 = P_Knee(1) + P_Knee(2)*cos(x*P_Knee(end)) + P_Knee(3)*sin(x*P_Knee(end)) + P_Knee(4)*cos(2*x*P_Knee(end)) + P_Knee(5)*sin(2*x*P_Knee(end)) + P_Knee(6)*cos(3*x*P_Knee(end)) + P_Knee(7)*sin(3*x*P_Knee(end)) + P_Knee(8)*cos(4*x*P_Knee(end)) + P_Knee(9)*sin(4*x*P_Knee(end)) + P_Knee(10)*cos(5*x*P_Knee(end)) + P_Knee(11)*sin(5*x*P_Knee(end));
    %         x2 = x + 0.5;
    %         diff = y1(end) - y3(1);
    %         y3 = y3 + diff;
    %         diff = y2(end) - y4(1);
    %         y4 = y4 + diff;
    %
    %         plot([x, x2],[y1, y3]);
    %         plot([x, x2],[y2, y4]);
    %
    %     end
    
end







