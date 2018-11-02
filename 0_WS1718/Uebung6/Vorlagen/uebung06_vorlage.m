% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Übung 5:
% Kontinuierliche und Diskrete Differentialoperatoren
%
% Inhalt:
% - Numerische Berechnung räumlicher Ableitungen
%

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('5. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 1 %%%
% %%%%%%%%%%%%%%%%%

% Aufgabe 1.a): Erstellung eines äquidistanten Gitters:
Nsteps = [ 5, 10, 20, 40, 80, 1000 ]; % verschiedene Gitterweiten ausprobieren
for k = 1:length(Nsteps)
    N     = Nsteps( k ); % Aktuellen Gitterparameter wählen
    Nvert = (N+1)*(N+1); % Nvertices: Anzahl der Knoten im Netz
    X     = zeros(Nvert, 1); % X-Koordinaten
    Y     = zeros(Nvert, 1); % Y-Koordinaten
    for i =             % i: Index in x-Richtung
        x = 
        for j =         % j: Index in y-Richtung
            y      = 
            idx    = VertexIndex(i,j,N);
            X(idx) = x; % Koordinaten
            Y(idx) = y; % setzen
        end
    end

    %%% Knoten als Kreise darstellen:
    plot(X, Y, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    %%% Bedeutet:
    %%% Plotte die Paare (X_i,Y_i) mit Symbol 'o' (s. Ü4-Folien)

    % Aufgabe 1.b):
    % Ableitung nach der x-Koordinate:
    DX = DiscreteDX( N ); % Diskreter Differentialoperator als Matrix
    % Bandstruktur von DX visualisieren:
    . . . . .

    % Ableitung nach der y-Koordinate:
    DY = DiscreteDY( N );
    % Bandstruktur von DY visualisieren:
    . . . . .
    
    % Aufgabe 1.c):
    % Der Differentialoperator DIFF:
    DIFF =  . . . . .
    % Bandstruktur von DIFF visualisieren:
    . . . . . 
    
        
    % Aufgabe 1.d):
    % LAPLACE OPERATOR:
    LAPLACE =  . . . . .
    % Bandstruktur von LAPLACE visualisieren:
    . . . . . 
    wait = waitforbuttonpress 
    
    % Aufgabe 1.e):
    % Vergleich zwischen diskreter und analytischer Lösung:
    f    = zeros(Nvert,1);
    dxf  = zeros(Nvert,1);
    dyf  = zeros(Nvert,1);

    % f1:
    f    = exp(- X .* Y);
    dxf  = -Y .* f;
    dyf  = -X .* f;
    lapf = ( Y .* Y + X .* X) .* f;
    divf = dxf + dyf;

    % f2:
%     f    = cos( 10*X .* Y);
%     dxf  = -  10 * Y .* sin( 10*X .* Y);
%     dyf  = -  10 * X .* sin( 10*X .* Y);
%     lapf = - 100 * (Y .* Y + X .* X ) .* f;
%     difff = dxf+dyf;

    % Fehlerberechnung:
    errdxold  = 1.;
    errdyold  = 1.;
    errdivold = 1.;
    if ( k > 1 )
        errdxold  =  errdx;
        errdyold  =  errdy;
        errdiffold = errdiff;
        errlapold = errlap;
    end
    errdx   = . . . . .
    errdy   = . . . . .
    errdiff = . . . . .
    errlap  = . . . . .
    disp(sprintf('h        = %22.14g',   1 / N ));
    disp(sprintf('ERR_dx   = %22.14g',   errdx ));
    disp(sprintf('ERR_dy   = %22.14g',   errdy ));
    disp(sprintf('ERR_diff = %22.14g', errdiff ));
    disp(sprintf('ERR_lap  = %22.14g',  errlap ));
    if( k > 1 ) % Nur für k>1 mit dem alten Fehler vergleichen
        disp(sprintf('ln (ERR_dx^new/ERR_dx^old)   / ln(h^new/h^old) = %22.14g', log10( errdx /   errdxold ) / log10( Nsteps(k-1)^2/ N^2 ) ));
        disp(sprintf('ln (ERR_dy^new/ERR_dy^old)   / ln(h^new/h^old) = %22.14g', log10( errdy /   errdyold ) / log10( Nsteps(k-1)^2/ N^2 ) ));
        disp(sprintf('ln (ERR_div^new/ERR_diff^old)/ ln(h^new/h^old) = %22.14g', log10(errdiff/ errdiffold ) / log10( Nsteps(k-1)^2/ N^2 ) ));
        disp(sprintf('ln (ERR_lap^new/ERR_alp^old) / ln(h^new/h^old) = %22.14g', log10( errlap/  errlapold ) / log10( Nsteps(k-1)^2/ N^2 ) ));
    end

end

% Nur letzte Lösung darstellen, da besonders "interessant" (besser: exakt)
clf;
hold all;
plot3(X,Y,. . . . .,  'o', 'Color','b');
plot3(X,Y,. . . . .*f, 'd', 'Color','r');
legend('exakt','numerisch');
