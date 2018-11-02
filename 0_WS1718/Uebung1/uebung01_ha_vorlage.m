% Rechnerunterstützte Mechanik I
% Wintersemester 2017/18
%
% - Schleifen- und Verzweigungsprogrammierung (IF ...)
% - Formatierte Bildschirmausgabe auf die Matlab-Konsole
% 

disp('Hausaufgabe zur 1. Übung im Fach Rechnerunterstützte Mechanik I (WS 2017/18)');

% Schleifen mit Verzweigung:
k       = 0.0; % Wert der aktuellen Zufallszahl
kmin    = 0.75;
counter = 0;
j       = 0;
Nnumber = 10;
disp(sprintf('So lange Zufallszahlen auswürfeln, bis %d-mal eine Zufallszahl > %f erreicht wird' , Nnumber, kmin));
while ( . . . . )
    k = 
    
    %%% Formatierte Ausgabe %%%
    % -> MATLAB Hilfe verwenden!
    % 18 Stellen insgesamt, Fließkommazahl
    disp(sprintf('X=%18f',k));
    % 22 Stellen insgesamt, wissenschaftliche Notation, 16 Nachkommastellen
    % disp(sprintf('X=%22.16g',k));

    % Verzweigungsprogrammierung: Zufallszahl größer als 
    if ( . . . . ) 
        . . . . % Anzahl der gefundenen Variablen erhöhen
    end
    
    counter = . . . .  % Counter für Anzahl Zufallszahlen inkrementieren
end

% %d : Dezimalzahl
disp(sprintf('Es wurden %d Versuche benötigt um %d Werte größer 0.75 zu "würfeln".', counter, nfound));

