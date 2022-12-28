function [xkplusone] = discreteDyn(xk, uk, params)
%discreteDyn State transition for linear discrete dynamics

xkplusone = params.A * xk + params.B * uk + params.d;
end