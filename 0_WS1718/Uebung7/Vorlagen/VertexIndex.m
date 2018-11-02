function i = VertexIndex(xi, yi, N)
% Gibt den Index eine Knotens zurück.
% ACHTUNG: Wertebereich: xi,yi = 1..N+1
    i = xi + (yi-1) * (N + 1);
end