% RechnerunterstÃ¼tzte Mechanik I
% Inhalt:
% - Iterativen LÃ¶ser selber programmieren (CG,Steepest Descent)
%
clc;clear all;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Gradientenverfahren   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B       = [3,2;2,6];
y       = [2;-8];
errmax  = 1.e-8;
normy   = norm(y,2);
err     = 1;
NMAX    = 100;
niter   = 0;                   % Anzahl der Iterationen
N       = size(B,1);
x       =     ....             % LÃ¶sungsvektor
f       = 1/2*x'*B*x-x'*y;     %Funktionswert der quadratischen Funktion für den Startwert
r       =    ....              %Residumm/Abstiegsrichtung in erster Iteration
PSIG(1,:)=[x;f];                % Abspeichern von aktueller Lösung und Funktionswert zur späteren Visualisierung
while(       )
    %Bestimmen Sie die Lösung x der aktuellen Iteration
    %Berechnen Sie für die Visualisierung den Funktionswert f der
    %quadratischen Funktion
    % Gradientenverfahren siehe Übung
 PSIG(niter+1,:)=[x;f];
end
disp(sprintf('Gradientenverfahen'))
disp(sprintf('Konvergenz nach %i Iterationen mit Residuum %e',niter,norm(r,2)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   CG-Verfahren       %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matrix einlesen 
% Verwenden Sie zuerst die gegebene 2x2 Matrix aus dem Gradientenverfahren
% prüfen Sie hiermt ob die Verfahen funktionieren
%Bzero = importdata('sparse_matrix2.mat');
%k     = 2;
%B=Bzero*Bzero';
niter    = . . . . % Anz. Iterationen, 
NMAX     = 150; % Max. Anz. Iterationen,
x        = . . . . % LÃ¶sungsvektor
r        = . . . . % Residuum
rold     = . . . . % Residuum der vorigen Iteration

d        = zeros(N,1);
p        = zeros(N,1);

err      = ;       %Bestimmung des relativen Fehlers ||Bx_j-y||/||y||
errmax   = 1.e-8;
normy    = max(1.e-10, norm(y,2));

RR       = r'*r;   % Spart Rechenzeit
RRold    = 0.;     % Spart Rechenzeit
f        = ....... %Funktionswert der quadratischen Funktion für Startwert
PSICG(1,:)=[x;f]; % Abspeichern von aktueller Lösung und Funktionswert zur späteren Visualisierung 
while( . . . . . . )

    %%% Neue Abstiegsrichtung (Fallunterscheidung fÃ¼r erstes Inkrement!)
    if ( . . . . )

    
    else
        
        
    end
    
    . . . .
    . . . .
    x      = . . . .
    rold   = r;
    r      = . . . .
    err    = . . . .
    
   %Ausrechnen des Funktionswert der quadtratischen Funktion
   f=
   %Abspeichern der aktuellen x-Werte und des Funktionswertes in einer
   %Matrix zur Visualisierung
    PSICG(niter+1,:)=[x;f];
    disp(sprintf('Iteration   %3d, relative residual   %16.10g,  ||Bx_j-y||/||y||= %16.10g',niter, err, norm(B*x-y,2)/normy));
    
end

%Oberflächenplot der quadratischen Form
x1=-5:0.05:5;
x2=x1;
for i=1:size(x1,2);
    for j=1:size(x2,2);
        xplot=[x1(i),x2(j)];
        f(i,j)=1/2*xplot*B*xplot'-xplot*y;
    end
end
figure
hold all;
grid on;
surf(x2,x1,f)
%Plot des Lösungsvektors für jede Iteration des jeweiligen Verfahren
%grün Gradient 
%rot  CG-Verfahren
plot3(PSIG(:,2),PSIG(:,1),PSIG(:,3),'o-g')
plot3(PSICG(:,2),PSICG(:,1),PSICG(:,3),'o--r')



