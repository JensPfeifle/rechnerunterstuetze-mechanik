% Rechnerunterstützte Mechanik I
% Wintersemester 2016/17
%
% Übung 3:
% Große Lineare Gleichungssystem
%
% Inhalt:
% - LU-Zerlegung (in Anlehnung an Dolittle-Algorithmus)

clc;
disp('3. Übung zur Vorlesung Rechnerunterstützte Mechanik I (WS 2016/17');

%%%%%%%%%%%%%%%%%
%%% AUFGABE 1 %%%
%%%%%%%%%%%%%%%%%

A = [ 
    ];

% Dimension von A:
N = 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1.a)    - Dolittle-Algorithmus (LUP-Zerlegung)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ LU, P ] = LUPFaktorisierung( A );

%%% Die Matrix LU in L und U zerlegen:
L  = eye(N);  U  = eye(N);
. . . . . . . .

%%% Residuum (in Matrixform) berechnen:
R = 
disp('Frobenius-Norm des Residuums:');
disp(sprintf('||LUP-A||_F :  %22.14g\n', . . . . ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.b)     - Gleichungslösung:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Zufällige rechte Seite:
y =
% Aufruf des selbstgeschriebenen Solvers mit Fehlerkorrektur:
x = SolveWithLUP(L, U, P, y);

% Residuum:
r = 

% Norm von r ausgeben:
disp('Residuum:');
disp(sprintf('||Ax-y||1   :  %22.14g', . . . . . ));
% DITO: ||.||2  ;   ||.||max

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.c)     - Test mit zufälligen Matrizen:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndim = [ 10, 20 . . . . . .

for ii=1:length(ndim)
    N = ndim(ii);
    A = . . . . .

    [LU,P] = . . . . .
    R      = . . . . .
    disp(sprintf('Dimension: %5d,  || A-LUP ||F = %22.14g', . . . . . ));
end

%%%%%%%%%%%%%%%%%
%%% AUFGABE 2 %%%
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.a)     - Test mit vorgegebener Matrizen:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix B einlesen:

B =  

% Struktur der Matrix ausgeben lassen:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.b)     - GAUSS-SEIDEL-Algorithmus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bzero = 
N     = 
% Symmetrisieren:
B     = 

% Diagonaldominanz herstellen:
k = 

% Ober und untere Dreicksmatrix initialisieren:
LplusD = 
U      = 


EPSILON = % Tolaranzparameter
n       = % Anzahl Iterationen
e       = % Relativen Fehler initialisieren
y       = % Rechte Seite des GLS
x       = % Lösungsvektor
Nmax    = % max. Anzahl Iterationen
Noutput = 5;
while(  )

    % Residuum berechnen:
    r = B*x - y;
    % Relativer Fehler:
    e = 

    % Ausgabe des Residuums alle Noutput Iterationen (zur Konvergenzkontrolle)
    if(rem(n,Noutput) == 0 )
        . . . . .
    end
end
% Ausgabe Anzahl Iterationen:




