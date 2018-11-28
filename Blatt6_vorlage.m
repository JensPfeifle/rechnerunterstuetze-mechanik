% Rechnerunterstuetzte Mechanik I
% Wintersemester 2018/19
%
% Uebung 5:
% Kontinuierliche und Diskrete Differentialoperatoren
%

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('6. Uebung im Fach Rechnerunterstuetzte Mechanik I (WS 2018/19)');

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 1 %%%
% %%%%%%%%%%%%%%%%%
disp('Aufgabe 1');

% Frequenz zur Berechnung von omega (für Randbedingungen)
freq    = 1;
omega   = 2*freq*acos(-1);
% Amplitude (für Randbedingungen)
amp  = 100;

% Funktion zur Definition der Randwerte
Randwert = @(x,y) amp * ( sin ( omega * x ) + sin ( omega * y )) + 293;
    
% Verschiedene Schrittweiten definieren:  
Nsteps = [ 5, 10, 20, 40 ]; % verschiedene Gitterweiten ausprobieren

for k = 1:length(Nsteps)
    N = Nsteps( k ); % Aktuellen Gitterparameter wählen

    N     = Nsteps( k ); % Aktuellen Gitterparameter wählen
    Nvert = (N+1)*(N+1); % Nvertices: Anzahl der Knoten im Netz
    X     = zeros(Nvert, 1); % X-Koordinaten
    Y     = zeros(Nvert, 1); % Y-Koordinaten
    % ToDo: Erstellung eines äquidistanten Gitters

    % Temperatur initialisieren
    theta = zeros(Nvert,1);

    % ToDo: Randwerte setzen
    
    % ToDo: Rechte Seite für unsere Gleichung setzen

    % Dimension der reduzierten Matrix
    Nvertmod= (N-1)*(N-1);
    
    % ToDo: Modifikation der Systemmatrix
    
    % Werte der Sparse Matrix auslesen
    [ i, j, val ] = find(LAPLACE);
    nval = size(val,1); % Anzahl der von Null verschiedenen Werte

    % Indizes für die modifizierte Matrix initialisieren
    imod = zeros(1, nval);
    jmod = zeros(1, nval);
    vmod = zeros(1, nval);
    
    % ToDo: [imod, jmod, vmod] aus [i,j,v] erstellen
    
    % Modifizierte Sparse Matrix erstellen
    Lmod = sparse( imod(1:nmod), jmod(1:nmod), vmod(1:nmod), Nvertmod, Nvertmod);
    
    % ToDo: Rechte Seite des Gleichungssystems anpassen
    
    % Konditionszahl der modifizierten Matrix
    disp(sprintf('Kondition der modifizierten Systemmatrix %22.14g\n', ...
                        condest(Lmod)));     % Konditionszahl

    % ToDo: Lösung des Gleichungssystems mittels des bicg - Lösers

    % Die Lösung liegt in der reduzierten Form vor
    % ToDo: Rücktransformation für theta
    
    % ToDo: Ausgabe der Lösung:
    % siehe auch surfc Hilfe
    % Eingabe: ToDo: Liste der x-Koordinaten des regulären Gitters MyX
    %          ToDo: Liste der y-Koordinaten des regulären Gitters MyY
    %          (Achtung: Nicht Punktepaare, sondern alle x
    %                    und alle y Werte)
    %          ToDo: Matrix der Temperaturwerte MyT

    surfc( MyX, MyY, MyT); 
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
    randwerte(ct) = Randwert(x,y);
    randX(ct)     = x;    randY(ct)     = y;    ct = ct + 1;
end
% rechts:
for j = 0:N-1
    x = 1;    y = j/N;
    randwerte(ct) = Randwert(x,y);
    randX(ct)     = x;    randY(ct)     = y;    ct = ct + 1;
end
% oben:
for i = N:-1:1
    x = i/N;    y = 1;
    randwerte(ct) = Randwert(x,y);
    randX(ct)     = x;    randY(ct)     = y;    ct = ct + 1;
end
% links:
for j = N:-1:1
    x = 0;    y = j/N;
    randwerte(ct) = Randwert(x,y);
    randX(ct)     = x;    randY(ct)     = y;    ct = ct + 1;
end
% ... und zurück zum Start:
randwerte(4*N+1) = randwerte(1);
randX(4*N+1) = randX(1);
randY(4*N+1) = randY(1);
plot3(randX, randY, randwerte, 'Color','m','LineWidth',2);
colorbar;
legend off;
