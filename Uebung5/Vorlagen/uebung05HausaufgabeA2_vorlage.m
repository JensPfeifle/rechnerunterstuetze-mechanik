% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Übung 4: Zusatzaufgabe
% Differentialoperatoren
% Anwendung auf gewöhnliche Differentialgleichungen
%

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('Hausaufgabe 2 zur 4. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

% %%%%%%%%%%%%%%%%%%%%%%%
% %%% HAUSAUFGABE 2 %%%
% %%%%%%%%%%%%%%%%%%%%%%%

%
% Implizites Lösungsverfahren
%

% Vektor mit den verschiedenen Zeitschrittanzahlen:
Ntsteps = [ 5, 20, 100, 500 ]; % Anzahl der Zeitschritte
rho     =   % Dichte
c       =   % spez. Wärmekapazität
doth    =   % dot{h}_0
k       =   % k
t0      =   % Startzeit
Tmax    =   % t_max
qab     =   % Abwärme
EPSTOL  =   % Toleranz für den nicht-linearen Löser (klein genug!)

hold all;        % Plot-Overlay
for i=1:length(Ntsteps)
    % Zeitintervall anpassen:
    N  = 
    dt = 
    T  = 
    
    % Theta initialisieren:
    theta    = 
    % Anfangswert:
    theta(1) = 
    
    for j=1:N
        % Residuum mit großem Wert initialisieren
        r = 1;
        % Expliziter Prädiktorschritt:
        dtheta     = dt* ( . . . . . );
        theta(j+1) = % ACHTUNG: dtheta wird in der while-Schleife addiert!
        niter      = % Anzahl Iterationen initialisieren
        % Maximal 10 Iterationsschritte (sonst: nächster Zeitschritt)
        while ( abs(r) > EPSTOL && niter < 10)
            % Theta korrigieren
            theta(j+1) = 
            % Neues Residuum berechnen:
            r          =     theta(j+1) - . . . . .

            % Korrekturterm neu berechnen:
            dfx        = . . . . . . .
            dtheta     = -r/(1 - dt*dfx);
            % Anzahl Iterationen inkrementieren:
            niter      = niter+1;

            % Zur Kontrolle des Residuums auskommentieren:
            % disp(sprintf('Residuum (Iteration %4d): %22.14g',niter,abs(r)));
            
        end
    end
    % Temperaturverlauf darstellen
    
end
% Legende und Titel hinzufügen
legend( . . . . . );
title( . . . . ,'FontSize',14,'FontWeight','bold')
