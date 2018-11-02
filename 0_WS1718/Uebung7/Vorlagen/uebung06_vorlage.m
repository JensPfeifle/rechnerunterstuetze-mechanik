% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Übung 6:
% Lösung der stationären Wärmegleichung mit Dirichletrandbedingungen
%

 . . . . .    % Konsole löschen
 . . . . .    % Alle Figures löschen & schließen
 . . . . .    % Speicher leeren

disp('6. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

% Frequenz zur Berechnung von omega (für Randbedingungen)
freq    = 2;
omega   = . . . . .
% Amplitude (für Randbedingungen)
amp  = 100;

% Verschiedene Schrittweiten definieren:
Nsteps = [ . . . . . ];
    
% FOR-Schleife über alle Diskretisierungsvarianten:
for k = 1:length(Nsteps)
    N     =         % Aktuelle Diskretisierung wählen
    Nvert =         % Anzahl der Knoten initialisieren
    % Speicherplatz reservieren:
    X     =         % x-Koordinaten ALLER Punkte
    Y     =         % dito y-Koordinaten
    % Werte von X und Y setzen
    for i = 1:(N+1)                 % i: y-Koordinate
        x = (i-1)/N;                % für festes i ist y konstant
        for j = 1:(N+1)             % j: x-Koordinate
            idx    = VertexIndex(i, j, N);
            y      = (j-1)/N;
            X(idx) = x;
            Y(idx) = y;
        end
    end

    % Funktion zur Definition der Randwerte
    % in Abhängigkeit von x und y definieren
    Randwert = @(x,y) amp * ( sin ( omega * x ) + sin ( omega * y )) + 293;

    % LAPLACE
    LAPLACE =  . . . . . .;

    % Temperatur initialisieren
    theta = zeros(Nvert,1);

    % Randwerte setzen:
    for i = 1:(N+1)
        for j = 1:(N+1)
            % Auf dem Rand von Null verschiedene Werte definieren:
            % || : Logisches ODER
            % && : Logisches UND
            if ( ( i == 1 ) || . . . . .  || ( j == N+1 ) )
                idx = VertexIndex(i, j, N);
                x   = X(idx);
                y   = Y(idx);
                theta(idx) = . . . . . . ;
            end
        end
    end

    % Rechte Seite für unsere Gleichung setzen:
    rhs   = . . . . . ;

    % Diskreten Operator ändern (wg. Randwerten):
    Nvertmod= (N-1)*(N-1);
    % Speicherplatz reservieren:
    rhsNew  = zeros( Nvertmod, 1 );

    % Modifikation der Systemmatrix
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ein schneller Algorithmus ist notwendig;
    % Deshalb wird hier die "sparse" Eigenschaft
    % der Matrix LAPLACE ausgenutzt:
    [ i, j, val ] = find(LAPLACE);
    nval = size(i,1); % Anzahl der von Null verschiedenen Werte

    % Indizes für die modifizierte Matrix anlegen
    % ACHTUNG: Später wird nicht der gesamt index verwendet, um
    %          die sparse-Matrix zu erzeugen
    imod = zeros(1, nval);
    jmod = zeros(1, nval);
    vmod = zeros(1, nval);
    nmod = 0;
    for ii=1:nval
        % Index xI1, ... reichen von 1..N+1
        % Zeilenindex in x- und y-Komponente unterteilen:
        [ xI1, yI1 ] = InverseVertexIndex( i(ii), N ); 
        [ xI2, yI2 ] = InverseVertexIndex( j(ii), N );
        % Prüfen, ob ein innerer Punkt vorliegt
        if (  ( xI1 > 1 ) && ( xI1 <= N ) && ...
                ............................
              ( yI2 > 1 ) && ( yI2 <= N )  )
            nmod = nmod+1;
            imod(nmod) = VertexIndex(xI1 - 1, yI1 - 1, N - 2);
            jmod(nmod) = VertexIndex(xI2 - 1, yI2 - 1, N - 2);
            vmod(nmod) = val(ii);
        end
    end
    Lmod = sparse( imod(1:nmod), jmod(1:nmod), vmod(1:nmod), Nvertmod, Nvertmod);
    % Rechte Seite des Gleichungssystems muss ebenfalls angepasst werden:
    for i = 2:N                                    % Schleife über innere Punkte
        for j = 2:N                                % Schleife über innere Punkte
            rowOld = VertexIndex(   i,   j,   N ); % Index des Punktes in
                                                   % der originalen Matrix
            rowNew = VertexIndex( i-1, j-1, N-2 ); % Index des Punktes in
                                                   % der neuen Matrix
            rhsNew( ......... ) = ..... ;       % Rechte Seite permutieren
        end
    end

    % Bandstruktur visualisieren
    . . . . 

    % Konditionszahl der modifizierten Matrix berechnen (kleiner = besser)
    disp(sprintf('Kondition der modifizierten Systemmatrix %22.14g\n', ...
                        condest(Lmod)));     % Konditionszahl

    % Lösung des Gleichungssystems:
    disp(sprintf('Löse das System mit Hilfe des Matlab Biconjugate Gradient Gleichungslöser.\n\n'));

    % Aufruf des bicg - Lösers:
    MySolution = bicg( . . . . . , 1e-8, 1000);

    % Die Lösung liegt in der reduzierten Form vor
    % (innere Punkte des Gitters)
    % -> Rücktransformation:
    for i = 2 : N
        for j = 2 : N
            rowOld = VertexIndex( .... );
            rowNew = VertexIndex( .... );
            theta( ... ) = ....;
        end
    end
    % Ausgabe der Lösung:
    % siehe auch surfc Hilfe
    % Eingabe: Liste der x-Koordinaten des regulären Gitters
    %          Liste der y-Koordinaten des regulären Gitters
    %          (Achtung: Nicht Punktepaare, sondern alle x
    %                    und alle y Werte)
    MyX = 0:1/N:1;
    MyY = 0:1/N:1;
    MyT = zeros(N+1, N+1); % Matrix mit den Temperaturwerten
    for i = 1:(N+1)
        for j = 1:(N+1)
            idx = VertexIndex( i, j, N );
            MyT( . . . . . ) = theta( . . . . );
        end
    end
    surfc( MyX, MyY, MyT); % OPTIONAL: 'EdgeColor', 'none'
    % zum Ausblenden des Gitters
    colormap jet
    
    k=waitforbuttonpress;
end


hold on;
legend('Temperatur \fontnames{times}\theta');


% Darstellung der Randwerte in einer anderen Farbe und als Linienzug:
randwerte = zeros(4 * N + 1 , 1);
randX     = zeros(4 * N + 1 , 1);
randY     = zeros(4 * N + 1 , 1);
ct        = 1;
% unten:
for i = 0:N-1
    x = i/N;    y = 0;
    randwerte(ct) = amp*(sin(omega*x) + sin(omega*y)) + 293;
    randX(ct)     = x;    randY(ct)     = y;    ct = ct + 1;
end
% rechts:
for j = 0:N-1
    x = 1;    y = j/N;
    randwerte(ct) = amp*(sin(omega*x) + sin(omega*y)) + 293;
    randX(ct)     = x;    randY(ct)     = y;    ct = ct + 1;
end
% oben:
for i = N:-1:1
    x = i/N;    y = 1;
    randwerte(ct) = amp*(sin(omega*x) + sin(omega*y)) + 293;
    randX(ct)     = x;    randY(ct)     = y;    ct = ct + 1;
end
% links:
for j = N:-1:1
    x = 0;    y = j/N;
    randwerte(ct) = amp*(sin(omega*x) + sin(omega*y)) + 293;
    randX(ct)     = x;    randY(ct)     = y;    ct = ct + 1;
end
% ... und zurück zum Start:
randwerte(4*N+1) = randwerte(1);
randX(4*N+1) = randX(1);
randY(4*N+1) = randY(1);
plot3(randX, randY, randwerte, 'Color','m','LineWidth',2);
colorbar;
legend off;
