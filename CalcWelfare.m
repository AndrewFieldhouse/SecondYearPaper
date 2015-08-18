function [W,CW, LTW] = CalcWelfare(M_W,gamma,sigma)
% Calculates welfare based on given allocation
W = zeros(402,1);
CW = zeros(402,1);
c = M_W(:,1);
l = 1 - M_W(:,2);
d = M_W(:,3);

% Toggle 1 for seperate public goods, 0 for 
%sw = 0;
sw =1;

wedge = 0.15;

if length(M_W(1,:)) == 4
    
    g = M_W(:,4);
    
    % Perfect Substitutes
    for i = 1:402        
        W(i) = ((d(i)^(i-1))*(((c(i)^gamma)*l(i)^(1-gamma))^(1-sigma)-1)/(1-sigma)+0.18*log(g(i)));
        if i == 1
            CW(i) = W(i);
        else
            CW(i) = CW(i-1) + W(i);
        end
    end
elseif length(M_W(1,:)) == 5
    
    gf = M_W(:,4);
    gs = M_W(:,5);
    
    % Perfect Compliments
    if sw == 0
        for i = 1:402
            W(i) = ((d(i)^(i-1))*(((c(i)^gamma)*l(i)^(1-gamma))^(1-sigma)-1)/(1-sigma)+0.18*log(2*min(gf(i), gs(i))));
            if i == 1
                CW(i) = W(i);
            else
                CW(i) = CW(i-1) + W(i);
            end
        end
        
    % Totally distinct public goods
    elseif sw == 1

        for i = 1:402
            W(i) = ((d(i)^(i-1))*(((c(i)^gamma)*l(i)^(1-gamma))^(1-sigma)-1)/(1-sigma)+0.09*log(gf(i))+0.09*(1-wedge)*log(gs(i)));
            if i == 1
                CW(i) = W(i);
            else
                CW(i) = CW(i-1) + W(i);
            end
        end
    end
end
LTW = sum(W(1:21));