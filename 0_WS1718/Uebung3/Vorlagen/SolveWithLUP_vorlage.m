function [ x ] = SolveWithLUP( L, U, P, y )

EPSILON =  % Relative Toleranz
NMAX    =  % Max. Iterationen
N       =  % Dimension von L/U/P

n =  % Anzahl Iterationen?
r =  % Residuum
x =  % Lösungsvektor
z =  % Zwischenvektor für Vor-/Rückwärtseinsetzen
v =  % Zwischenvektor für Permutation
e =  % Relativer Fehler

while( . . . . . . )
    % Lz = y lösen nach z:

    % Uv = z lösen nach v:
    
    % x aktualisieren:
    
    % Zähler inkrementieren:
    n=n+1;
    
    
    % Residuum aktualisieren
    
    % Relativen Fehler berechnen:

% Testweise Residuum ausgeben lassen:
%    disp(sprintf('Residuum: (iteration %d)',n));
%    disp(sprintf('||Ax-y||2 :  %22.14g', e));
end
