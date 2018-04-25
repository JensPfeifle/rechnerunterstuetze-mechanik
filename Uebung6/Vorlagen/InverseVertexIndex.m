function [ xi, yi ] = InverseVertexIndex( k, N )
xi = rem( k-1, N+1 ) + 1;
yi = (k-xi)/(N+1) + 1;
end