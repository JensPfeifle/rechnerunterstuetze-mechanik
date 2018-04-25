function plotEmodFibo(C, phiN) 
% Input: C: Steifigkeit in Abaqus Notation
%	 phiN: Anzahl der jeweiligen Unterteilungen der beiden Winkel zur Beschreibung der Raumrichtung (Achtung: Gesamtanzahl = phiN^2)
%
% Output: Plot des E-Modul-Körpers

% Transformation der Steifigkeit von Abaqus Notation in normierte Voigt-Notation
Dsqrtt  = diag([ 1, 1, 1,  sqrt(2),  sqrt(2),  sqrt(2) ]);
C2      = Dsqrtt*C*Dsqrtt;

% Unterteilung aller Raumrichtungen in phiN^2 Punkte ( phi1 in [0,2*pi], phi2 in [0,pi] )
%phi1=linspace(0,2*pi,phiN);
%phi2=linspace(0,pi,phiN);

% Unterteilung aller Raumrichtungen mittels Fibonacci Verteilung
phi=(1+sqrt(5))/2;
i=(-(phiN-1):2:(phiN-1))';
theta =2*pi*i/phi;
sphi=i/phiN;
cphi=sqrt((phiN+i).*(phiN-i))/phiN;

n = zeros ( phiN, 3 );
n(1:phiN,1) = cphi .* sin ( theta );
n(1:phiN,2) = cphi .* cos ( theta );
n(1:phiN,3) = sphi;

S = inv(C2);

%DTV = zeros(3,phiN^2);
DTV = zeros(3,phiN);
DTVi = 1;

hold all;
%for i=1:size(phi2,2)
%    for j=1:size(phi1,2)
for j=1:phiN
	% Aufstellen der Normalenrichtung für einachsigen Zug
        %nv = ([cos(phi1(j)).*sin(phi2(i)); sin(phi1(j)).*sin(phi2(i)); cos(phi2(i))]);
        nv = n(j,:);

	% Aufstellen von n dyadisch n
        NdN=([nv(1)^2; nv(2)^2; nv(3)^2; sqrt(2)*nv(2)*nv(3); sqrt(2)*nv(1)*nv(3); sqrt(2)*nv(1)*nv(2)]);
        
	% Berechnung des Betrags des E-Moduls
        ET=1/trace(NdN*(S*NdN)');
        
	% E-Modul in Normalenrichtung
        DTV(:,DTVi) = ([ET*nv(1);ET*nv(2);ET*nv(3)]);
        
        DTVi = DTVi+1;
%    end
end

% 3d plot über alle Winkel
plot3(DTV(1,:),DTV(2,:),DTV(3,:),'.');

axis equal
