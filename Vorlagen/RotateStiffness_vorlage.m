function [Crot , RR ] = RotateStiffness( . . . , . . . , . . . , . . . ) 

% Basiswechsel:
Dsqrth  = diag([ 1, 1, 1,  1/sqrt(2),  1/sqrt(2),  1/sqrt(2) ]);
Dsqrtt  = diag([ 1, 1, 1,  sqrt(2),  sqrt(2),  sqrt(2) ]);							

% Aufstellen der Rotationsmatrix
s1   = sin(phi1);   s2   = sin(PHI);    s3   = sin(phi2);
c1   = cos(phi1);   c2   = cos(PHI);    c3   = cos(phi2);

% 3x3 Rotationsmatrix
Q = zeros(3);
Q = [ . . . ] 								

% Orthonormalbasis für normierte Voigt-Notation
B = zeros(6,3,3);
fac  =  1/sqrt(2);
B(1, . . . , . . . ) = . . . ;
B(2, . . . , . . . ) = . . . ;
B(3, . . . , . . . ) = . . . ;
B(4, . . . , . . . ) = . . . ;
B(5, . . . , . . . ) = . . . ;
B(6, . . . , . . . ) = . . . ;

% Orthonormalbasis im gedrehten System
Brot = zeros(6,3,3);
for i=1:6
   Brot(i,:,:) = . . . ;
end

RR = zeros(6,6);
for i=1:6
    for j=1:6
        RR(i,j) = . . . ;
    end
end

% Rotierten Steifigkeitstensor berechnen
Rabq = . . . ;
Crot = . . . ;