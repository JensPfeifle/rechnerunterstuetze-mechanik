% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% Übung 4:
% Differentialoperatoren
%

clc;       % Konsole löschen
close all; % Alle Figures löschen & schließen
clear;     % Speicher leeren

disp('4. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 1 %%%
% %%%%%%%%%%%%%%%%%

% 1.a)
N     =  
t0    =  
tN    = 
T     = 

% Definition einer feinen Unterteilung zur Darstellung
% einer Referenzlösung
Nfine = 
Tfine = 


% 1.b)
% Funktionsdefinition:
% s. f1.m, f2.m, f3.m
%
% Eingabe für f1, f2, f3: Zeitpunkt(e) (s. Aufgabenstellung)
% Ausgabe von f1, f2, f3: Funktionswert(e)

% 1.c) Differentialoperatoren:
% s. ForwardEuler.m, BackwardEuler.m, MidpointRule.m
%
% EINGABEPARAMETER:
% i.   Indexmenge, z.B. [ 0 1 5 7 12 ]
% ii.  Vektor der diskreten Zeitpunkte T
% iii. Zeiger auf die Funktion, die verwendet werden soll
%      (s. Beispiel unten). Diese Funktion soll als Eingabeparameter
%      Zeitpunkt(e) t akzeptieren und einen Vektor von Funktionswerten
%      zurückgeben.
% RÜCKGABEWERTE:
% Vektor der durch das entsprechende Verfahren approximierten Ableitung
%
% Zu beachten:
% - Die Indexmenge kann Werte von 0 bis N enthalten.
% - Berechnen Sie in Ihrem Differenzenverfahren die Zeitschrittlänge für
%   jeden Wert. Damit ist das Programm auch für nicht gleichgroße
%   Zeitschritte einsetzbar.
% - Beachten Sie den Indexshift: t_0 <-> T(1), t_1 <-> T(2), etc.


% 1.d)
% Indexmenge initialisieren:
idx = 0:1:N; % Alternativ z.B. [ 0 1 5 7 12 ]

% Figure löschen
clf;
hold all; % Von jetzt an alles übereinander plotten


plot(T, BackwardEuler(idx, T, @f1), 'Color', 'red');
% Testen für
% ... ForwardEuler ...
% ... MidpointRule ...
% ... Exakte Lösung ...

% Auf Tastendruck/Mausklick warten:
k = waitforbuttonpress;

clf;
% wie oben für Funktion 2 und 3

% Alle figures schließen und Speicher leeren
close all;
clear;

% %%%%%%%%%%%%%%%%%
% %%% AUFGABE 2 %%%
% %%%%%%%%%%%%%%%%%

% 2.a) Systemmatrix für Vorwärts-Euler:
N  = 100;
t0 = 0;
tN = 1;

deltaT =
AFW    = sparse(zeros(N+1,N+1)); % dünnbesetzte Matrix erstellen
for i=1:N
    ...
    ...
end
% Eine Zeile fehlt noch ...


% 2.b) Systemmatrix für Rückwärts-Euler:
ABW    = . . . .
for i=2:N+1
    ...
    ...
end
% Eine Zeile fehlt noch ...

% 2.c) Systemmatrix für Mittelpunktsregel:
AMP    = . . . .
for(i=2:N)
end
% Erste und letzte Zeile modifizieren..

% %%%%%%%%%%%%%%%%%%%%%%
figure; % Neue Figure erstellen:
gcf;    % GetCurrentFigure
clf;    % ClearFigure

% Bandstruktur darstellen 
....
% Titel setzen:
title('Darstellung der dünnen Bandmatrix AMP','FontSize',14,'FontWeight','bold');
% Achsenbeschriftung:
set(get(gca,'XLabel'),'String','Spalte','FontSize',12);
set(get(gca,'YLabel'),'String','Zeile','FontSize',12);
% weitere Infos: Matlab Hilfe 'text'

% Funktionswerte als Vektor berechnen:
func1  = . . . . ;

% Diskrete Werte der Ableitungen über Matrix-Vektor-Produkt:
dfunc1 = . . . . ;


figure; % Neue Figure erstellen:
gcf;    % GetCurrentFigure
clf;    % ClearFigure
plot( - ForwardEuler -  'Color', 'red');
plot(Tfine, df1(Tfine)); % Darstellung der Referenzlösung
set(get(gca,'XLabel'),'String','{\fontname{times}\itt} [s]','FontSize',12);
set(get(gca,'YLabel'),'String','\fontname{times}{\itf}''_1(\itt)','FontSize',12);
hold all;
% BW-Euler: BLAU; Mittelpunktregel: GRÜN

% Legende, Titel, Achsenbeschriftung und exakte Lösung hinzufügen

func2 = f2(T)';
% .... und func3 analog!


