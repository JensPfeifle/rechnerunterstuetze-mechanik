function value=VectorNorm(x, p)
% Funktionsbeschreibung:
% Berechnet ||x||p = ( sum_i=1^n |x|^p ) ^ (1/p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Länge des Vektors
n     =   
% Rückgabewert initialisieren
value =   

% Die Sonderfälle p=1 und p=2 sind numerisch effizienter zu berechnen, als
% der allgemeine Fall für reelles p.
if ( p == 1 )
    % 1-Norm

elseif ( p == Inf )
    % Maximums-Norm

elseif ( p == 2 )
    % 2-Norm (Euklid-Norm)
    
    % Verwenden Sie vektorisierte Kommandos (KEINE FOR-Schleifen!)

else
    % p-Norm

end