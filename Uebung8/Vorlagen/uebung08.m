% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Übung 8:
% Elastizitätstheorie
%
% Inhalt:
% - Steifigkeitstensor
% - Hooke'sches Gesetz
% - von Mises Vergleichsspannung

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('8. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 1 %%%
% %%%%%%%%%%%%%%%%%

% 1.a) s. StressStrain.m
% 1.b) s. SigmaMises.m

% 1.c) s. CISO.m, CCUBIC.m, CTISO.m

% 1.d)
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

% 1.d) s. SpectralDecomposition.m
[ P, lam , nlam ] =  SpectralDecomposition( C );
CTest = zeros(6,6);
for i=1:nlam
   CTest = CTest + lam(i)*squeeze(P(i,:,:)); 
end
norm(CTest-C,'fro')

plotEmod(C, 1000);
