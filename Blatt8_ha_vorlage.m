% Rechnerunterstützte Mechanik I
% Wintersemester 2018/19
%
% Uebung 8:
% Numerische Integration
%
% Inhalt:
% - Newton-Cotes-Formeln
% - Gauss-Quadratur (eindimensional)

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('Testat zur 8. Uebung zur Vorlesung Rechnerunterstuetzte Mechanik I');

omega = 13.6753;
PI    = acos(-1);
Tmin  = 0;
Tmax  = 0.5*PI;
funk  = @(x)  sin(omega.*x);

% 1.b) TODO: (handschriftlich bestimmen!) Exakter Wert des Integrals:
Iexakt = . . . 

%%% Gauss-Christoffel-Quadratur:
disp('');
disp('Gauss-Quadratur');
%%% Wie viele Koeffizienten sollen berechnet werden?
maxorder = 10;
N        = maxorder;
Legendre = zeros(N + 1, N + 1);
Nst      = zeros(N + 1, N + 1); %%% Dim=Ordnung + 1

%%% Berechnung der Legendre-Polynome 1...maxorder+1:
%%% (die NST des N+1-ten Polynoms liefern die Stuetzstellen
%%% fuer das Verfahren der Ordnung N)
for i=1:N+1
    % TODO: Koeffizienten alpha ausrechnen mittels IntegratePolynomial
    alpha = zeros(i,1);
     
    %%% TODO: Ausgehend von p=t^(i-1) das Legendre-Polynom aufstellen

    %%% TODO: Polynom normieren

end

%%% Wenig aussagekraeftige Ausgabewerte...
disp('Koeffizienten der Legendre-Polynome:');
Legendre

%%% TODO: Nullstellenberechnung mittels roots Befehl (liefert Stuetzstellen)
%%% Speichern in Matrix tau
tau = zeros(N,N); % Speicher reservieren

disp('');
disp('Stützstellen für die Gauss-Quadratur:');
tau
save('gauss_tau.mat','tau');

% TODO: Berechnung der Lagrange-Gewichte:
% Speichern in Matrix lambda
lambda=zeros(N,N);
for i=1:N % i: Anzahl der Stuetzstellen des Verfahrens
    % T: Vektor der Stützstellen
    T(1:i) = tau(i,1:i);
    for j=1:i
        p    = zeros(i,1);
        % Starte mit p(t) = 1, d.h. p(1) = 1;
        p(1) = 1;
        
        % TODO: Durchlaufe alle Stuetzstellen zur Berechnung des Polynoms

        % TODO: Berechnung des Lagrangegewicht mittels IntegratePolynomial
    end
end
disp('');
disp('Gewichte:');
lambda
save('gauss_lambda.mat','lambda');

load('errorNC.mat')
% Auswertung des Integrals (ein Teilintervall):
for i = 1:N
    T   = zeros(i,1);
    for j = 1:i
        T(j) = Tmax*tau(i,j);
    end
    f   = funk(T);
    val = Tmax*lambda(i,1:i)*f;
    disp(sprintf('Approximation mit %2d Stützstellen: %22.14g, Fehler: %22.14g, Fehler_NC: %22.14g', i, val, Iexakt-val, errorNC(i)));
end

disp('');
disp('Berechnung mit mehreren Teilintervallen:');

% Grafische Darstellung der Funktion sowie der Stützstellen
n = 100;
x = Tmin:(Tmax-Tmin)/n:Tmax;
figure;
plot(x,funk(x))
hold on

% Auswertung des Integrals (mehrere Teilintervall):
Nint     = 10; % Anzahl der Teilintervalle
integral = 0.; % Initialisierung des Integrals
order    = 4;  % Anzahl der Stützstellen PRO TEILINTERVALL
h        = (Tmax-Tmin)/Nint; % Schrittweite
T        = zeros(order,1); % Vektor der Stuetzstellen im aktuellen Teilintervall
for k=0:Nint-1
    % TODO: Definition der Intervallgrenzen a, b

    % TODO: Bestimmung des Funktionswerts an der Stuetzstelle

    % TODO: Bestimmung des Gewichts
    
    % TODO: Integral erhoehen
    
    plot(T,f,'*');
end
disp(sprintf('Approximation mit %2d Teilintervallen mit je %2d Stützstellen: %22.14g, Fehler: %22.14g', Nint, order, integral, Iexakt-integral));

line(xlim,[0 0],'Color','r')
