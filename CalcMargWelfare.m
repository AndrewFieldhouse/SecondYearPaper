function [margW] = CalcMargWelfare(M_W,gamma,sigma)
% Calculates welfare based on given allocation
%W = zeros(402,1);
%CW = zeros(402,1);
c = M_W(1,1);
l = 1 - M_W(1,2);
%d = M_W(:,3);
%g = M_W(:,4);

margW = (((c^gamma)*(l)^(1-gamma))^(-sigma))*gamma*(c^(gamma-1))*(l^(1-gamma));