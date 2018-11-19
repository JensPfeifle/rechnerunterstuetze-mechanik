% Rechnerunterstuetzte Mechanik I
% Wintersemester 2018/19
%
% Uebung 4:
% Differentialoperatoren
%

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('4. Uebung im Fach Rechnerunterstuetzte Mechanik I (WS 2018/19)');

%%%%%%%%%%%%%%%%%
%%% AUFGABE 1 %%%
%%%%%%%%%%%%%%%%%
disp('Aufgabe 1');

% Aufgabe 1a) 
N     = 500; 
Nfine = 1000;
t0    = 0;
tN    = 1;

% ToDo: Definition de dimensionslosen Vektoren T und Tfine

% Aufgabe 1b) 
% ToDo: Schreiben Sie die Funktionen
% ForwardEuler(i, T, f)
% BackwardEuler(i, T, f)
% MidpointRule(i, T, f)
% welche die numerischen Ableitungen einer Funktion f berechnen.
% Hierbei soll ein Vektor der diskreten Zeitpunkte T uebergeben werden.
% Die Rueckgabe ist der Vektor der durch das entsprechende Verfahren
% approximierten Ableitung.

% Aufgabe 1c) 
% ToDo: Gemeinsame Darstellung mit der exakten Loesung

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 2 %%%
% %%%%%%%%%%%%%%%%%
disp('Aufgabe 2');

% Aufgabe 2.a) 
AFW    = sparse(zeros(N+1,N+1)); % duennbesetzte Matrix initialisieren
% ToDo: Systemmatrix für Vorwaerts-Euler:

% Aufgabe 2.b) 
AFB    = sparse(zeros(N+1,N+1)); % duennbesetzte Matrix initialisieren
% ToDo: Systemmatrix für Rueckwaerts-Euler:

% Aufgabe 2.c)
AMP    = sparse(zeros(N+1,N+1)); % duennbesetzte Matrix initialisieren
% ToDo: Systemmatrix für Mittelpunktsregel:

% Aufgabe 2.d)
% ToDo: Funktionswerte als Vektor berechnen

% ToDo: Diskrete Werte der Ableitungen über Matrix-Vektor-Produkt

% ToDo: Geinsame Darstellungen der Approximationen mit der exakten Loesung

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 3 %%%
% %%%%%%%%%%%%%%%%%
disp('Aufgabe 3');

% Vektor mit den verschiedenen Zeitschrittanzahlen:
Ntsteps = [ 1, 5, 20, 100, 500 ];
rho     = 1000;     %[kg/m^3]
c       = 1000;     %[J/kg*K]
doth    = 1.7e7;    %[J/s*m^3]
k       = 561;      %[K]
t0      = 0;
Tmax    = 100;
qab     = 7.5e5;  %[J/s*m^3]

% Figure löschen und hold aktivieren
figure;
clf;
hold on;
for i=1:length(Ntsteps)
    N  = Ntsteps(i);    % Anzahl Zeitschritte für das i-te Element aus Ntsteps
    dt = (Tmax-t0)/N;
    T  = t0:dt:Tmax;
    
    % Theta initialisieren
    theta    = zeros(N+1,1);
    theta(1) = 293;     %[K]
    
    % TODO: Berechnung von theta
    
    plot(T, theta); % akutelles Ergebnis in das bestehende Schaubild zeichnen
end
% Legende und Titel hinzufügen
% Titel setzen:
title('Wärmeentwicklung im Reaktor','FontSize',14,'FontWeight','bold');
legend('1 Zeitschritt', '5 Zeitschritte', '20 Zeitschritte', '100 Zeitschritte', '500 Zeitschritte');
% Achsenbeschriftung:
set(get(gca,'XLabel'),'String','{\fontname{times}\itt} [s]','FontSize',12);
set(get(gca,'YLabel'),'String','\theta [K]','FontSize',12);


% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 4 %%%
% %%%%%%%%%%%%%%%%%
disp('Aufgabe 4');

EPSTOL  = 1e-10;  % Toleranz für den nicht-linearen Löser (klein genug!)

% Figure löschen und hold aktivieren
figure;
clf;
hold on;
for i=1:length(Ntsteps)
    % Zeitintervall anpassen:
    N  = Ntsteps(i);
    dt = (Tmax-t0)/N;
    T  = t0:dt:Tmax;
    
    % Theta initialisieren
    theta    = zeros(N+1,1);
    theta(1) = 293;     %[K]
    
    % TODO: Berechnung von theta
    for j=1:N
        niter = 0; % Anzahl Iterationen initialisieren
        
        % ToDo: Expliziter Prädiktorschritt:
        
        % ToDo: Residuum berechnen
        
        % Maximal 10 Iterationsschritte (sonst: nächster Zeitschritt)
        while ( abs(r) > EPSTOL && niter < 10)
            % ToDo: Theta korrigieren, Residuum berechnen, Iteration inkrementieren
            
            % Zur Kontrolle des Residuums
            disp(sprintf('Residuum (Iteration %4d): %22.14g',niter,abs(r)));
        end
    end
    % Temperaturverlauf darstellen
    plot(T, theta); % akutelles Ergebnis in das bestehende Schaubild zeichnen
end
% Legende und Titel hinzufügen
title('Wärmeentwicklung im Reaktor','FontSize',14,'FontWeight','bold');
legend('1 Zeitschritt', '5 Zeitschritte', '20 Zeitschritte', '100 Zeitschritte', '500 Zeitschritte');
