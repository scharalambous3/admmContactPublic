function [xkplusone] = discreteDyn(xk, uk, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% xkplusone = zeros(size(xk, 1), 1);

xkplusone = params.A * xk + params.B * uk + params.d;
end