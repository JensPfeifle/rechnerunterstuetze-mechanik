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

disp('Hausaufgabe 1 zur 4. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

% %%%%%%%%%%%%%%%%%%%%%%%
% %%% HAUSAUFGABE 1 %%%
% %%%%%%%%%%%%%%%%%%%%%%%

% 1.a) / 1.b)

Ntsteps = [ . . . . . ];
rho     = 
c       = 
doth    = 
k       = 
t0      = 
Tmax    = 
qab     = 0;

% Figure löschen und hold aktivieren
for i=1: . . . .
    N  = % Anzahl Zeitschritte für das i-te Element aus Ntsteps

    dt = 
    T  = 
    
    % Theta initialisieren
    theta    = 
    theta(1) = 
    
    for j=1:N
        theta(j+1) = % Iterationsvorschrift
    end
    plot(T, theta); % akutelles Ergebnis in das bestehende Schaubild zeichnen
end
% Legende und Titel hinzufügen
legend('1 Zeitschritt', '5 Zeitschritte', '20 Zeitschritte', '100 Zeitschritte', '500 Zeitschritte');
title('Wärmeentwicklung im Reaktor','FontSize',14,'FontWeight','bold')
% Achsenbeschriftung:
set(get(gca,'XLabel'),'String','{\fontname{times}\itt} [s]','FontSize',12);
set(get(gca,'YLabel'),'String','\theta [K]','FontSize',12);



