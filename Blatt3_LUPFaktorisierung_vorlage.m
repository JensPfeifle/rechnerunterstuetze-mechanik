function [ LU, P ] = LUPFaktorisierung( A )
% Berechnet die LU-Zerlegung der Matrix A
% Es gilt:
% L = LU(i,j) ( i > j ); 1 ( i == j); 0 (sonst)
% U = LU(i,j) ( j >= i ); 0 (sonst)

% Speicherplatz reservieren:
N  = size(A, 1);
LU = zeros(N, N);

% Eingabematrix kopieren
LU = A;

% Timer initialisieren
t0 = cputime;

% ToDo: LU-Zerlegung berechnen:
P = speye(N);
npivot = 0; % Anzahl der Pivotisierungen
for i = 1:N  % i: Zeile
    for j = 1:N  % j: Spalte
        . . . . .
    end

    % ToDo: Spalten-Pivotsuche:
    % Maximales Element in der aktuellen Zeile
    . . . . . 

    % ToDo: Permutation der Spalten durchfuehren:
    if( . . . . . ) % Nur, falls Permutation notwendig
        % ToDo: Vertausche LU(:,i) mit LU(:,n) und P(i,:) mit P(n,:):
        . . . . .
    end  
end
Rechenzeit = cputime-t0;
disp(sprintf('Rechenzeit: %22.8f', Rechenzeit));
disp(sprintf('Spaltenpermutationen: %d (von %d Spalten)', npivot, N));
