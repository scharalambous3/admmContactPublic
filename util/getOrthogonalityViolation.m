function [viol] = getOrthogonalityViolation(Z_k, params)
%getOrthogonalityViolation Returns violation of orthogonality constraint of
%complementarity constraint

N = size(Z_k, 2) + 1;
viol = 0;
for i = 1:(N - 1)
    %TODO
    viol = viol + norm((params.Aorth * Z_k(:, i) + params.aorth)' * (params.Borth * Z_k(:, i) + params.borth));
end
end