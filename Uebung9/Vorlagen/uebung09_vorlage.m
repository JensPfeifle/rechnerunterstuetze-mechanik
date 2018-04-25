% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Übung 9:
% Numerische Integration
%
% Inhalt:
% - Newton-Cotes-Formeln
% - Gauss-Quadratur (eindimensional)
% - 
%

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('9. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 1 %%%
% %%%%%%%%%%%%%%%%%

% Newton Cotes - Formeln
disp('');
disp('Aufgabe 1 - Newton-Cotes-Quadratur');
% 1.a): Anfängliche Definitionen
omega = 2 ; 
PI    = acos(-1);     % Funktioniert auch in C++ oder Fortran am besten so!
Tmin  = . . . . . ;
Tmax  = . . . . . ;   % Intervallbreite festlegen
% Definition der Funktion aus 1.a)
% Beispiel zum Syntax: funk = @(x) x.*x; % Entspricht: (funk(x))_i = x_i^2
funk  = @(x) . . . . . ;

% 1.b): (handschriftlich bestimmen!) Exakter Wert des Integrals:
Iexakt = . . . . . ;

% 1.c): s. 'IntegratePolynomial.m'
% TESTEN:
p = [ 0, 1, 0, 0, -1 ]; % Beispielpolynom p(t) = t - t^4;
% Test:
IntegratePolynomial( p, Tmin, Tmax );

% Wie viele Koeffizienten sollen berechnet werden?
Ncoeffs = 10;
coeffs  = zeros(Ncoeffs, Ncoeffs+1);
for i=1:Ncoeffs
    Nc     = i;
    T      = 0:Tmax/Nc:Tmax;
    for j=1:Nc+1
        p    = zeros(Nc+1,1);    % Genügend Speicherplatz für das Polynom reservieren
        p(1) = . . . . . ;                % Start: p(t) = 1
        for k=1:Nc+1             
            if(j ~= k)
                . . . . . ;
                % Schritt 1: p nach rechts shiften:
                for l=Nc:-1:1
                    . . . . . ;
                end
                % Ersten Koeffizienten löschen:
                . . . . . ;
                % Schritt 2: Addieren von t_k*p_old:
                . . . . . ;
                % Schritt 3: Multiplizieren mit 1/(t_j-t_k):
                . . . . . ;
            end
        end % Insgesamt Nc-1 Aufrufe der if-Schleife
        % Koeffizient speichern:
        coeffs(i, j) = IntegratePolynomial(p, Tmin, Tmax);
        save('coeffs.mat', 'coeffs');
    end
end

% Koeffizienten einlesen und ausgeben:
coeffs = 0;
load('coeffs.mat')
disp('Koeffizientenmatrix');
coeffs

% Exakten Wert ausgeben:
disp(sprintf('Wert des exakten Integrals: %22.14g',Iexakt));

% Numerische Berechnung:
for i=1:Ncoeffs
    T= . . . . . . % Zeitintervall entsprechend teilen
    f= . . . . . . % Vektor der Funktionswerte an den Stützstellen
    val = . . . . . .  % Wert über einfaches Vektor-Matrix-Produkt
    % Anzahl Stützstellen, numerisches Integral und numerischen Fehler ausgeben:
    disp(sprintf('Approximation mit %2d Stützstellen: %22.14g, Fehler: %22.14g', . . . . . . . ));
end

% Grafische Darstellung der Funktion sowie der Stützstellen
n = 100;
x = Tmin:(Tmax-Tmin)/n:Tmax;
figure;
plot(x,funk(x))
line(xlim,[0 0],'Color','r') 
hold on
plot(T,f,'*');
