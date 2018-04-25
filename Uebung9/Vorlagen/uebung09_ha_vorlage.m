% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Übung 9:
% Numerische Integration
%
% Inhalt:
% - Gauss-Quadratur
%

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('Testat zur 9. Übung zur Vorlesung Rechnerunterstützte Mechanik I');

omega = 13.6753;
PI    = acos(-1);
Tmin  = 0;
Tmax  = 0.5*PI;
funk  = @(x)  sin(omega.*x);

%%% Exakter Wert des Integrals:
Iexakt = . . . . . ;

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
    % Koeffizienten alpha ausrechnen:
    alpha = zeros(i,1);
  
    %%% Berechnung der Koeffizienten:
    for j=1:i-1
        q = zeros(2*i,1); %%% Ausreichend Speicher reservieren
        %%% Shift right: (q entpricht: Legendre_j * t^(i-1) )
        q(i:i+j-1) = Legendre(j,1:j);
        %%% Berechne t^(i-1)*Legendre_j
        alpha(j)   = . . . . . ;
    end
    
    %%% Jetzt muss ausgehend von p=t^(i-1) das Legendre-Polynom
    %%% aufgestellt werden
    Legendre(i,i) = 1;
    for j = 1:i-1
        for k = 1:j
            Legendre(i, k) = . . . . . ;
        end
    end
    %%% Das Polynom muss nun noch normiert werden:
    NormFaktor = . . . . . ;
    Legendre(i,1:i) = . . . . . .;
end

%%% Wenig aussagekraeftige Ausgabewerte...
disp('Koeffizienten der Legendre-Polynome:');
Legendre

%%% Nullstellenberechnung: (liefert Stuetzstellen)
tau = zeros(N,N); % Speicher reservieren
for i = 1:N
    p = zeros(i+1,1);
    %%% ACHTUNG: Die Zeilen aus der Matrix
    %%% Legendre haben ein anderes Format, als fuer
    %%% den Matlab internen roots-Befehl vorgesehen
    %%% -> Umspeichern:
    for j = 1:i+1
        p(j) = . . . . . ;
    end
    %%% Nullstellenberechnung:
    r = roots( p );
    %%% Nullstellen in Matrix tau ablegen:
    tau(i,1:i) = . . . . ;
end

disp('');
disp('Stützstellen für die Gauss-Quadratur:');
tau
save( . . . . );

% Berechnung der Lagrange-Gewichte:
lambda=zeros(N,N);
for i=1:N % i: Anzahl der Stuetzstellen des Verfahrens
    % T: Vektor der Stützstellen
    T(1:i) = . . . .;
    for j=1:i
        p    = zeros(i,1);
        % Starte mit p(t) = 1, d.h. p(1) = 1;
        p(1) = 1;
        
        % Durchlaufe alle Stuetzstellen
        for k = 1:i
            % Bei der Berechnung des j-ten Polynoms muessen alle
            % Stuetzstellen ausser der j-ten beruecksichtigt werden
            if ( . . . . )
                pold = p; % Aktuelles Polynom zwischenspeichern
                % Koeffizienten von p um eins nach rechts
                % verschieben (entspricht Multiplikation mit t):
                for l=i:-1:2
                    p(l) = . . . .
                end
                % Der erste Koeffizient muss NULL sein:
                p(1) = 0;
                % p_neu = ( t*p_alt - tau_k*p_alt )/ ( tau_j - tau_k )
                p = . . . . . ;
            end
        end
        % Polynom konstruiert. Berechne Integral -> Lagrangegewicht
        lambda(i,j) = . . . . . ;
    end
end
disp('');
disp('Gewichte:');
lambda
save( . . . . );

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
Nint     = 5; % Anzahl der Teilintervalle
integral = 0.;
order    = 4;  % Anzahl der Stützstellen PRO TEILINTERVALL
h        = . . . . . ; % Schrittweite
T        = zeros(order,1); % Vektor der Stuetzstellen im aktuellen Teilintervall
for k=0:Nint-1
    a = . . .;  b = . . .; % Intervallgrenzen

    for j = 1:order
        T(j) = a + h*tau(order,j);
    end
    f   = . . . .
    val = . . . .

    integral = integral + val;
    
    plot(T,f,'*');
end
disp(sprintf('Approximation mit %2d Teilintervallen mit je %2d Stützstellen: %22.14g, Fehler: %22.14g', Nint, order, integral, Iexakt-integral));

line(xlim,[0 0],'Color','r')
