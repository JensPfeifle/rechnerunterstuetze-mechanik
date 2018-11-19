% Rechnerunterstuetzte Mechanik I
% Wintersemester 2018/19
%
% Uebung 3:
% Grosse Lineare Gleichungssystem
%
% Inhalt:
% - LUP Faktorisierung
% - CG Verfahren
%

clc;
disp('3. Uebung zur Vorlesung Rechnerunterstuetzte Mechanik I (WS 2018/19)');

%%%%%%%%%%%%%%%%%
%%% AUFGABE 1 %%%
%%%%%%%%%%%%%%%%%
disp('Aufgabe 1');

A = [ 
         0,   5,   0,   4;
        -1,  -6,   0,  -2;
         9,   3,  -2,   1;
        -4,   7,  -1,   3
    ];

% ToDo: Aufgabe 1a) LUP-Faktorisierung
% Schreiben Sie zunaechste eine Funktion
% LUPFaktorisierung( A )
% welche LUP Faktorisierung einer Matrix A beliebiger Dimension durchfuehrt.
% Die Matrizen L und U sollen in einer gemeinsamen Matrix zurueckgegeben werden. 
% Zusaetzlich soll die Permutationsmatrix P berechnet und zurueckgegeben werden.


% Aufgabe 1b) Loesung mittels LUP Faktorisierung
% ToDo: Schreiben Sie eine Funktion
% SolveWithLUP( L, U, P, y )
% welche die Loesung des Gleichungsystems A x = y mittels der berechneten LUP Faktorisierung durchfuehrt.
% Zerlegen Sie hierzu zunaechst die obigen Matrizen in L, U und P und
% pruefen Sie die Zerlegung mittels Residuum.
% Pruefen Sie ebenfalls den Fehler durch den gefundenen Loesungsvektors x

N = size(A,1);
y = rand(N,1);

%%%%%%%%%%%%%%%%%
%%% AUFGABE 2 %%%
%%%%%%%%%%%%%%%%%
disp('Aufgabe 2');

% Matrix B einlesen:
B = importdata('sparse_matrix.dat');
% Struktur der Matrix ausgeben lassen:
spy(B);
N = size(B,1);

A=B*B';
y = rand(N,1);

% ToDo: Loesen Sie das Gleichungssystem Ax=y mittels des CG Verfahrens, so
% dass das maximale Residuum errmax innerhalb NMAX Iterationsschritten
% erreicht wird
NMAX     = 150;
errmax   = 1.e-8;

% ToDo: Algorithmus initialisieren

% ToDo: While-Schleife mit Restriktion an Residuum und Anzahl der
% Schleifendurchgaenge

% ToDo: Ausgabe wie viele Iterationen benoetigt wurden und welches Residuum
% erreicht wurde

% ToDo: Vergleich der gefundenen Loesung und Iterationen mit Matlabs bicg Loeser