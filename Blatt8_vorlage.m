% Rechnerunterstuetzte Mechanik I
% Wintersemester 2018/19
%
% Uebung 8:
% Numerische Integration
%
% Inhalt:
% - Newton-Cotes-Formeln
% - Gauss-Quadratur (eindimensional)


clc;       % Konsole loeschen
close all; % Alle Figures loeschen & schließen
clear;     % Speicher leeren

disp('8. Uebung im Fach Rechnerunterstuetzte Mechanik I (WS 2018/19)');

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 1 %%%
% %%%%%%%%%%%%%%%%%

% Newton Cotes - Formeln
disp('');
disp('Aufgabe 1 - Newton-Cotes-Quadratur');
% 1.a): Anfaengliche Definitionen
omega = 2 ; 
PI    = acos(-1);     % Funktioniert auch in C++ oder Fortran am besten so!
Tmin  = 0;
Tmax  = PI/2;   % Intervallbreite festlegen

% 1.a) TODO: Definition der Funktion
% Beispiel zum Syntax: funk = @(x) x.*x; % Entspricht: (funk(x))_i = x_i^2

% 1.b) TODO: (handschriftlich bestimmen!) Exakter Wert des Integrals:
Iexakt = . . . 

% 1.c) TODO: Schreiben Sie die Funktion IntegratePolynomial, welche fuer
% den Koeffizientenvektor das Integral des zugehoerigen Polynoms auf dem
% Intervall a bis b berechnet. 
% Eingangsparameter: Polynom p, untere Intervallgrenze a, obere Intervallgrenze b

% Beispielpolynom p(t) = t - t^4;
p = [ 0, 1, 0, 0, -1 ]; 
% TODO: Test des Integrals

% Anzahl der maximal zu berechenden Koeffizienten. Es werden nicht nur die
% Koeffizienten maximaler Ordnung gespeichert, um spaeter den Vergleich
% der Integrationen mit unterschiedlicher Koeffizientenzahl zu
% ermoeglichen.
Ncoeffs = 8;
coeffs  = zeros(Ncoeffs, Ncoeffs+1);
for i=1:Ncoeffs
    % Jeweilige Anzahl der Koeffizienten
    Nc = i;
    % Unterteilung des Intervalls
    T = Tmin:(Tmax-Tmin)/Nc:Tmax;
    
    % Schleife ueber einzelne Koeffizieten
    for j=1:Nc+1
    % TODO: Berechnung des p Vektors
    
    % TODO: Berechnung des Geiwchts mittels der IntegratePolynomial
    % Funktion. Speichern in coeffs Matrix.
    
    end
    
    % Koeffizientenmatrix speichern (fuer Testat und spaetere Uebungen)
    save('coeffs.mat', 'coeffs');
end

% Koeffizienten einlesen (unnoetig, jedoch Test fuer obige Speicherung) und ausgeben:
coeffs = 0;
load('coeffs.mat')
disp('Koeffizientenmatrix');
coeffs

% Exakten Wert ausgeben:
disp(sprintf('Wert des exakten Integrals: %22.14g', Iexakt));

% TODO: Numerische Berechnung:
for i=1:Ncoeffs
   
    % Anzahl Stuetzstellen, numerisches Integral und numerischen Fehler ausgeben:
    disp(sprintf('Approximation mit %2d Stuetzstellen: %22.14g, Fehler: %22.14g', . . . . . . . ));
end

% Grafische Darstellung der Funktion sowie der Stuetzstellen
n = 100;
x = Tmin:(Tmax-Tmin)/n:Tmax;
figure;
plot(x,funk(x))
line(xlim,[0 0],'Color','r') 
hold on
plot(T,f,'*');
