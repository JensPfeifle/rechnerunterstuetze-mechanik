function df = MidpointRule(i, T, f)
% Hinweise in uebung4_vorlage.m beachten!
N  = % Größe des Vektors T
Ni = % Anzahl der übergebenen Indizes, an denen df ausgewertet werden soll
df = % Rückgabevektor
for j=1:Ni % alle Indizes durchlaufen:
    idx = i(j);
    if (  )
        % "Normal"-fall
    elseif (  ) % Sonderfall, in dem ForwardEuler versagt
                % Entsprechendes anderes Verfahren aufrufen
    elseif (  ) % Anderer Sonderfall

    else % Gar nicht zulässige Indizes (z.B. -1 etc.)
        df(j) = Inf; % Error!
    end
end
% Midpoint rule
