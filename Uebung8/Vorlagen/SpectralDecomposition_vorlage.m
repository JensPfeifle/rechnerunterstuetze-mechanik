function [ P, lambda, nlambda ] = SpectralDecomposition( C )

% Basiswechsel:
Dsqrth  = diag([ 1, 1, 1,  1/sqrt(2),  1/sqrt(2),  1/sqrt(2) ]);
Dsqrtt  = diag([ 1, 1, 1,  sqrt(2),  sqrt(2),  sqrt(2) ]);

C2      = . . . . ; %%% C in der normierten Voigt-Notation

%%% Berechne die Eigenwerte bezüglich der originalen Basis:
[ U, V ] = eig( . . . . );

%%% Speicherplatz reservieren:
P       = zeros(6,6,6);  %%% Projektoren
lambda  = zeros(6,1);    %%% Eigenwerte

nlambda = 0;           % Anzahl verschiedener Eigenwerte initialisieren
u       = zeros(6,1);  % Eigenvektor (temporaerer Vektor)

for i = 1:6
    found = 0; %%% Index, so dass  V(i,i) = lambda(found)

    %%% Falls i = 1 -> nlambda = 1, lambda(nlambda) = V(1,1)
    if ( nlambda == 0 ) %%% Ersten Eigenwert auf jeden Fall übernehmen
        nlambda         = 1;
        lambda(1)       = V(. . . . , . . . . );
        found           = 1;
    else
        %%% Falls schon mind. ein Eigenwert bearbeitet wurde:
        %%% auf doppelte EW prüfen!
        for j=1:nlambda
            %%% Falls 'numerisch gleich' Index speichern
            if ( abs( . . . ) < 1.e-8*norm(lambda,Inf))
                found = . . . ;
            end
        end
        if ( found == 0 ) %%% Neuer Eigenwert gefunden
            nlambda         = . . . . ;
            lambda(nlambda) = V( . . . , . . . ); 
            found           = . . . . ;
        end
    end

    %%% u: Eigenvektor zu lambda(i)
    %%% !!! ACHTUNG:  Die u_j zu _vielfachen Eigenwerten_ sind
    %%% !!!           i.A. NICHT orthogonal -> Orthogonalisieren:
    u(:)   = U(:,i);
    A      = zeros(6,6);
    A(:,:) = P(found,:,:); %%% Projektor des EW (temporär)
    v      = u - A*u; %%% v = u - P_found u
    %%% Eigenprojektor um v times v inkrementieren
    %%% Die Division sorgt für die notwendige NORMIERUNG (||v|| = 1)
    A      = A + v*v'/(v'*v); 
    %%% Neuen Projektor abspeichern:
    P(found,:,:)=A(:,:);
end

% Zurückrechnen von normierter Voigt-Notation nach ABQ-Notation
for i=1:nlambda
    P(i,:,:) =  . . . *squeeze(P(i,:,:))* . . . ;
end


