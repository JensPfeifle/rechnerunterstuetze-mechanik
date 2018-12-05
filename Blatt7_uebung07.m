% Rechnerunterstuetzte Mechanik I
% Wintersemester 2018/19
%
% Uebung 7:
% Elastizitaetstheorie
%
% Inhalt:
% - Steifigkeitstensor
% - Hooke'sches Gesetz
% - von Mises Vergleichsspannung

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('7. Uebung im Fach Rechnerunterstuetzte Mechanik I (WS 2018/19)');

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 1 %%%
% %%%%%%%%%%%%%%%%%
disp('Aufgabe 1');

% Aufgabe 1.a) StressStrain.m

% Aufgabe 1.b) SigmaMises.m

K = 152564;
G = 26923;
C = CISO( K, G);

C1111 = 168000;
C1122 = 121000;
C2323 = 75000;
%C = CCUBIC( C1111, C1122, C2323 );

C1111 = 5000;
C3333 = 2000;
C2323 = 3500;
C1122 = 1500;
C1133 = 1000;
%C = CTISO( C1111, C3333, C2323, C1122, C1133 );

% 1.c) SpectralDecomposition.m

[ P, lam , nlam ] =  SpectralDecomposition( C );
CTest = zeros(6,6);
for i=1:nlam
   CTest = CTest + lam(i)*squeeze(P(i,:,:)); 
end
norm(CTest-C,'fro')

plotEmodFibo(C, 2000);

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 2 %%%
% %%%%%%%%%%%%%%%%%
disp('Aufgabe 2');

phi1 = 0;
PHI  = pi/4;
phi2 = 0;

Crot = RotateStiffness( C, phi1, PHI ,phi2 );

hold on
plotEmodFibo(C, 2000);
plotEmodFibo(Crot, 2000);