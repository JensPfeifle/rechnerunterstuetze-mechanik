function plotEmod(C, phiN) 
% Input: C: Steifigkeit in Abaqus Notation
%	 phiN: Anzahl der jeweiligen Unterteilungen der beiden Winkel zur Beschreibung der Raumrichtung (Achtung: Gesamtanzahl = phiN^2)
%
% Output: Plot des E-Modul-Körpers

% Transformation der Steifigkeit von Abaqus Notation in normierte Voigt-Notation
Dsqrtt  = diag([ 1, 1, 1,  sqrt(2),  sqrt(2),  sqrt(2) ]);
C2      = Dsqrtt*C*Dsqrtt;

% Unterteilung aller Raumrichtungen in phiN^2 Punkte ( phi1 in [0,2*pi], phi2 in [0,pi] )
phi1=linspace(0,2*pi,phiN);
phi2=linspace(0,pi,phiN);

S = inv(C2);

DTV = zeros(3,phiN^2);
DTVi = 1;

hold all;
for i=1:size(phi2,2)
    for j=1:size(phi1,2)
	% Aufstellen der Normalenrichtung für einachsigen Zug
        nv = ([cos(phi1(j)).*sin(phi2(i)); sin(phi1(j)).*sin(phi2(i)); cos(phi2(i))]);
	% Aufstellen von n dyadisch n
        NdN=([nv(1)^2; nv(2)^2; nv(3)^2; sqrt(2)*nv(2)*nv(3); sqrt(2)*nv(1)*nv(3); sqrt(2)*nv(1)*nv(2)]);
        
	% Berechnung des Betrags des E-Moduls
        ET=1/trace(NdN*(S*NdN)');
        
	% E-Modul in Normalenrichtung
        DTV(:,DTVi) = ([ET*nv(1);ET*nv(2);ET*nv(3)]);
        
        DTVi = DTVi+1;
    end
end

% 3d plot über alle Winkel
plot3(DTV(1,:),DTV(2,:),DTV(3,:),'.');
axis equal