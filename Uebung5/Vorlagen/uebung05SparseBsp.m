% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Zusatzmaterial zur 4. Übung:
% Verwendung von sparse-Matrizen in Matlab
%

% Konsole löschen:
clc;
disp('4. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');
disp('Beispiele zur Verwendung von sparse-Matrizen');

% N: Dimension der Matrix
N = 1000;

% Einfache NxN sparse-Matrize anlegen (alle Einträge 0):
A = sparse( [], [], [], N, N);

% Einträge manipulieren:
% !!! ACHTUNG: Im Allgemeinen nicht so auf die Matrix zugreifen !!!
A(1,3) = 1; A(20,27) = 5;
for i=20:500
    A(400,i) = 1;
end
% Struktur visualisieren:
subplot(1,3,1);
spy(A);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix mit den Einträgen von s.o. direkt erzeugen:

% 1.) idxi : Zeilennummern für alle Einträge
% 2.) idxj : Spaltennummern für alle Einträge
% 3.) val  : Wert der Komponente A_(idxi,idxj)

idxi = [ 1, 20, 400*ones(1,481) ];
idxj = [ 3, 27, 20:500 ];
val  = [ 1, 5, ones(1,481) ];
B = sparse( idxi ,idxj, val, N, N );

% Probe:
norm(A-B,'inf');

subplot(1,3,2);
spy(B)

% Matrix A aus Speicher löschen:
% clear A;
clear A,B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Einfache Matrix erstellen, die
% aufeinanderfolgende Einträge eines Vektors addiert:
% (A x)_i = x_i + x_(i+1)  (bzw. für i=N: x_N)

% 1.) idxi : Zeilennummern für alle Einträge
% 2.) idxj : Spaltennummern für alle Einträge
% 3.) val  : Wert der Komponente A_(idxi,idxj)

N=20; % So sieht man die Struktur besser
idxi = [ 1 : N, 1 : N-1 ];
idxj = [ 1 : N, 2 : N   ];
val  = ones(1,length(idxi));
A = sparse( idxi ,idxj, val, N, N );

subplot(1,3,3);
spy(A);
