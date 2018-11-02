% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Hausaufgabe 8 - MusterlösungÜbung
%
% Rotation eines Steifigkeitstensors in ABAQUS Notation
%

clc;        % Konsole löschen
clf;        % Alle Figures löschen & schließen
clear;      % Speicher leeren
close all;

disp('Hausaufgabe zur 8. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

phi1 = 0;
PHI  = pi/4;
phi2 = 0;

%C = CISO( 152564, 26923);
%C = CCUBIC( 168000, 121000, 75000 );
%C = CTISO( 5000, 1000, 2000, 1000, 1000 );

C = CCUBIC( 168.4e3, 124.4e3, 75.39e3 );

C
Crot = RotateStiffness( C, phi1, PHI ,phi2 );
Crot

hold on
plotEmodFibo(C, 2000);
plotEmodFibo(Crot, 2000);
