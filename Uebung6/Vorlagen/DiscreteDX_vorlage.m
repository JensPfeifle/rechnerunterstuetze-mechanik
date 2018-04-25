function DX = DiscreteDX( N )
Nvert  = (N+1)*(N+1);
h      = 1./N;
iidx   = zeros(2*Nvert, 1); % Speicherplatz für alle Werte reservieren
jidx   = zeros(2*Nvert, 1); % Insgesamt besitzt die Matrix genau zwei
values = zeros(2*Nvert, 1); % Einträge je Zeile (also je Knoten)

% Zunächst: iidx0, jidx0, values0 bestimmen, so dass
%           die sparse Matrix (iidx0, jidx0, values0, Nvert, Nvert)
%           die Ableitung d / d_x_i in der ERSTEN ZEILE berechnet
%
% Vorgegeben:    Erste und letzte Zeile
% Zu Ergänzen:   Alle Zwischenwerte
iidx0  = . . . .
jidx0  = . . . .
values0= 1./h * . . . .

% Jetzt die Matrix Blockweise kopieren:
% Die Vektoren iidx0, jidx0, values0 müssen in
% die Vektoren iidx,jidx,values an der entsprechenden Stelle
% eingesetzt werden (Index-Shift um N+1 pro Zeile).
% Hinweise aus der Übung beachten!
for j=0:N
    iidx( j*2*(N+1)+1  : ... ) = iidx0 + ... *ones(1, ... );
    jidx( j*2*(N+1)+1  : ... ) = jidx0 + ... *ones(1, ... );
    values(j*2*(N+1)+1 : ... ) = values0;
end

DX = sparse(iidx,jidx,values,Nvert,Nvert);
