%% Automatically generates plots of n random kicks

function [Time, Hip, Knee] = Kick_Generator_Math(Starting, Ending, Minimum, Length)


load('P.mat');
Per_Mult = 1;
Per = 10;
isvalid = 0;
while isvalid == 0
    Per_Mult = Per_Mult * 1.03;
    i = randi([1,length(P)]);
%     Dy1 = {};
%     Dy2 = {};
    P_Hip_Original = P{i,1};
    P_Knee_Original = P{i,2};
    y1 = [];
    y2 = [];
    X = [];
%     Dy1_ii = {};
%     Dy2_ii = {};
    for ii = 1:length(P{i,1})
        P_Hip = P_Hip_Original{ii};
        P_Knee = P_Knee_Original{ii};
        
        x = linspace(0,1/length(P{i,1}));
        y1_ii = P_Hip(1)*sin(P_Hip(2)*x+P_Hip(3)) + P_Hip(4)*sin(P_Hip(5)*x+P_Hip(6)) + P_Hip(7)*sin(P_Hip(8)*x+P_Hip(9)) + P_Hip(10)*sin(P_Hip(11)*x+P_Hip(12)) + P_Hip(13)*sin(P_Hip(14)*x+P_Hip(15)) + P_Hip(16)*sin(P_Hip(17)*x+P_Hip(18)) + P_Hip(19)*sin(P_Hip(20)*x+P_Hip(21)) + P_Hip(22)*sin(P_Hip(23)*x+P_Hip(24));
        y2_ii = P_Knee(1)*sin(P_Knee(2)*x+P_Knee(3)) + P_Knee(4)*sin(P_Knee(5)*x+P_Knee(6)) + P_Knee(7)*sin(P_Knee(8)*x+P_Knee(9)) + P_Knee(10)*sin(P_Knee(11)*x+P_Knee(12)) + P_Knee(13)*sin(P_Knee(14)*x+P_Knee(15)) + P_Knee(16)*sin(P_Knee(17)*x+P_Knee(18)) + P_Knee(19)*sin(P_Knee(20)*x+P_Knee(21)) + P_Knee(22)*sin(P_Knee(23)*x+P_Knee(24));
%         dy1_ii = gradient(y1_ii) ./ gradient(x);
%         dy2_ii = gradient(y2_ii) ./ gradient(x);
%         dy1_ii = [dy1_ii(1), dy1_ii(end)];
%         dy2_ii = [dy2_ii(1), dy2_ii(end)];
        
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
%         Dy1_ii = [Dy1_ii, dy1_ii];
%         Dy2_ii = [Dy2_ii, dy2_ii];
        x = x + x(end)*(ii-1);
        X = [X, x];
        
    end
%     Dy1 = [Dy1; Dy1_ii];
%     Dy2 = [Dy2; Dy2_ii];
    
    
    canpass = 0;
    while canpass == 0
        y1 = [];
        y2 = [];
        X = [];
%         Dy1_ii = {};
%         Dy2_ii = {};
        Per = Per*Per_Mult;
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
%                 dy1_ii = gradient(y1_ii) ./ gradient(x);
%                 dy2_ii = gradient(y2_ii) ./ gradient(x);
%                 dy1_ii = [dy1_ii(1), dy1_ii(end)];
%                 dy2_ii = [dy2_ii(1), dy2_ii(end)];
                if ii == 1
                    y1_ii_end = y1_ii(1);
                    y2_ii_end = y2_ii(1);
                end
                y1_offset = y1_ii(1) - y1_ii_end;
                y2_offset = y2_ii(1) - y2_ii_end;
                y1_ii = y1_ii - y1_offset;
                y2_ii = y2_ii - y2_offset;
                
                
                if max(y1_ii) > 170 || max(y2_ii) > 170 || min(y1_ii) < 60 || min(y2_ii) < 60
                    pass = 0;
                    Per = Per .* 0.95;
                else
                    pass = 1;
                    Per = 5;
                end
                
                
                if pass == 1
                    y1_ii_end = y1_ii(end);
                    y2_ii_end = y2_ii(end);
                    y1 = [y1, y1_ii];
                    y2 = [y2, y2_ii];
%                     Dy1_ii = [Dy1_ii, dy1_ii];
%                     Dy2_ii = [Dy2_ii, dy2_ii];
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
        
        maxdiff = 30;
        
        if y1(1) < Starting(1) - maxdiff || y1(1) > Starting(1) + maxdiff
        elseif y2(1) < Starting(2) - maxdiff || y2(1) > Starting(2) + maxdiff
        elseif y1(end) < Ending(1) - maxdiff || y1(end) > Ending(1) + maxdiff
        elseif y2(end) < Ending(2) - maxdiff || y2(end) > Ending(2) + maxdiff
        elseif min(y1) < Minimum(1) - maxdiff || min(y1) > Minimum(1) + maxdiff
        elseif min(y2) < Minimum(2) - maxdiff || min(y2) > Minimum(2) + maxdiff
        elseif max(y1) > max([Starting, Ending]) + (2*maxdiff) || max(y2) > max([Starting, Ending]) + (2*maxdiff)
        else
            isvalid = 1;
        end
        
        
    end
    
    %         figure(i);
    %         hold on;
    %         plot(X, y1);
    %         plot(X, y2);
    Time = X;
    Hip = y1;
    Knee = y2;
    
    
    
    
    
    
end

Time = Time .* Length;



end