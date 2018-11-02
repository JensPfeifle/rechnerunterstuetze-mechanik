function [ LU, P ] = LUPFaktorisierung( A )
% Berechnet die LU-Zerlegung der Matrix A
% Es gilt:
% L = LU(i,j) ( i > j ); 1 ( i == j); 0 (sonst)
% U = LU(i,j) ( j >= i ); 0 (sonst)

% Speicherplatz reservieren:
N  = size(A, 1);

% Eingabematrix kopieren
LU = . . . .

% Timer initialisieren
t0 = cputime;


% LU-Zerlegung berechnen:
% Permutationsmatrix als (sparse) Identitätsmatrix initialisieren:
P = speye(N,N);

% Anzahl der Spalten-Pivot-Vertauschungen
npivot = 0;

for i = 1:N  % i: Zeile
    for j = 1:N  % j: Spalte
        if ( j < i )    % Berechnung von L

        else            % Berechnung von U

        end
    end

    % Spalten-Pivotsuche:
    % Maximales Element in der aktuellen Zeile von U
    maxel = norm(LU(i,i:N),Inf);

    % n: Index der Spalte mit dem größten Element:
    n = i;
    valmax = abs(LU(i,i));
    % Spaltenpermutation überhaupt nötig?
    if( valmax < 0.25*maxel )
        % Suche Index des größten Elements
        for ii = (i+1):N
            tmp = abs( . . . . . . );
            if ( tmp > valmax )
                % Index (n) und valmax aktualisieren:
                . . . . .
                . . . . .
            end
        end
        % Spaltenvertauschungs-Zähler inkrementieren
        npivot = . . . .
    end

    % Permutation der Spalten durchführen:
    if( i ~= n ) % (aber nur falls notwendig)
        % Vertausche LU(:,i) mit LU(:,n) und P(i,:) mit P(n,:):
        tmpcol  = zeros(N,1); %tempor�rer Vektor zum Abspeichern der aktuellen Spalte von LU
        tmpcol  = LU(:,i);
        LU(:,i) =. . .              % Spaltenvertauschung der LU-Matrix
        LU(:,n) =. . . 
        tmprow  = zeros(1,N);  % Zeilenvertauschung der P-Matrix
        tmprow  = P( . . . );
        P( . . . )  = P( . . . );
        P(n,:)  = tmprow;
    end
end

% CPU-Zeit messen und ausgeben:
Rechenzeit = cputime-t0;
disp(sprintf('Rechenzeit.............: %22.8f', Rechenzeit));
disp(sprintf('Spaltenpermutationen...: %6d (von %d Spalten)', npivot, N));
