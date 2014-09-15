%% Momento Invariante de Hu
% Es un promedio ponderado de la intensidad de cada pixel
% Propiedades: Area, centroide, orientacion

clear all
clc
%% Image reading 
Im = imread('lena.png');

Im=(rgb2gray(Im));
Im = imresize(Im,0.5);
% Im = imrotate(Im,30);

[row col] = size(Im)

Im = double(Im);

%% Raw moments
% Mij = Sumx Sumy [x^i * y^j * I(x,y)]
% I(x,y) -> intensidad del pixel
% In some cases considering image as a probability density funcion...
% ... Divide Mij by: Sumx Sumy I(x,y)

%Raw Moment:

    % function valor = raw_moment(i,j,fila,columna,matriz)
    % acum = 0;
    % for x=1:fila
    %     for y=1:columna
    %         acum = acum + (x^i)*(y^j)*matriz(x,y);
    %     end
    % end
    % valor = acum;

M00 = raw_moment(0,0,row,col,Im);   % i, j, fila, columna, imagen
M01 = raw_moment(0,1,row,col,Im);
M10 = raw_moment(1,0,row,col,Im);
M11 = raw_moment(1,1,row,col,Im);
M20 = raw_moment(2,0,row,col,Im);
M21 = raw_moment(2,1,row,col,Im);
M22 = raw_moment(2,2,row,col,Im);
M02 = raw_moment(0,2,row,col,Im);
M21 = raw_moment(2,1,row,col,Im);
M12 = raw_moment(1,2,row,col,Im);
M22 = raw_moment(2,2,row,col,Im);
M30 = raw_moment(3,0,row,col,Im);
M31 = raw_moment(3,1,row,col,Im);
M03 = raw_moment(0,3,row,col,Im);
M13 = raw_moment(1,3,row,col,Im);
    
%% Centroid's components
% Area for binary images or sum of grey level for greytone
% Sum: M00
% Centroid (x_prom, y_prom)
x_prom = M10/M00;  
y_prom = M01/M00;

%% Central moments: Translational Invariant
% mpq = Sumx Sumy [(x-x_prom)^p * (y-y_prom)^q * I(x,y)]
m00 = M00;
m01 = 0;
m10 = 0;
m11 = M11 - y_prom * M10;
m20 = M20 - x_prom * M10;
m02 = M02 - y_prom * M01;
m21 = M21 - 2*x_prom*M11 - y_prom*M20 + 2*x_prom^2*M01;
m12 = M12 - 2*y_prom*M11-x_prom*M02+2*y_prom^2*M10;
m30 = M30 - 3*x_prom*M20+2*x_prom^2*M10;
m03 = M03 - 3*y_prom*M02+2*y_prom^2*M01;

%% Orientation: Covariance Matrix
m20_p = m20 / m00;
m02_p = m02 / m00;
m11_p = m11 / m00;

% The covariance matrix of the image I(x,y) is:
cov = [m20_p m11_p;m11_p m02_p];
% The orientation can be extracted from the angle of the eigenvector...
% ... associated with the largest eigenvalue
angle_orient = 0.5 * atan(2*m11_p / (m20_p - m02_p));

%% Scale invariant moments
% The mpq moments can be constructed to be invariant to both translation...
%...  and change in scale. Divide the central moment by the scaled.

% n(i,j) = m(i,j)/(m(0,0)^(1+(i+j)/2))
n00 = m00/m00; 
n01 = m01/m00^(1.5);
n10 = m10/m00^(1.5); 
n11 = m11/m00^(2);
n20 = m20/m00^(2);
n02 = m02/m00^(2);
n21 = m21/m00^(2.5);
n12 = m12/m00^(2.5); 
n30 = m30/m00^(2.5);
n03 = m03/m00^(2.5);

%% Rotation invariant moments
% This is invariant to translation, scale and roatation: 
% Hu invariant moment

i1 = n20+n02;    % Moment of inercia arround the image's centroid
i2 = (n20-n02)^2+4*(n11)^2;
i3 = (n30-3*n12)^2+(3*n21-n03)^2;
i4 = (n30+n12)^2+(n21+n03)^2;
i5 = (n30-3*n12)*(n30+n12)*((n30+n12)^2 -3*(n21+n03)^2)+(3*n21-n03)*(n21+n03)*(3*(n30+n12)^2-(n21+n03)^2);
i6 = (n20-n02)*((n30+n12)^2-(n21+n03)^2)+4*n11*(n30+n12)*(n21+n03);
i7 = (3*n21-n03)*(n30+n12)*((n30+n12)^2-3*(n21+n03)^2)-(n30-3*n12)*(n21+n03)*(3*(n30+n12)^2-(n21+n03)^2);
i8 = n11*((n30+n12)^2-(n03+n21)^2)-(n20-n02)*(n30+n12)*(n03+n21);
% Sugieren que el i3 no es my úitl y agregan un octavo momento.

%% log transform

matriz_hu =  zeros(1,7,'double')

matriz_hu(1) = -sign(i1)*log10(abs(i1));
matriz_hu(2) = -sign(i2)*log10(abs(i2));
matriz_hu(3) = -sign(i3)*log10(abs(i3));
matriz_hu(4) = -sign(i4)*log10(abs(i4));
matriz_hu(5) = -sign(i5)*log10(abs(i5));
matriz_hu(6) = -sign(i6)*log10(abs(i6));
matriz_hu(7) = -sign(i7)*log10(abs(i7));
matriz_hu(8) = -sign(i8)*log10(abs(i8));
% Si la imagen es rotada, la comparación entre ésta y la original...
% ... producen una diferencia de signo en algunos coeficientes...
% ... Para la comparación idealmente ocmparar las magnitudes de tal modo...
% ... de aislar el error.
plot(matriz_hu,'yellow')
hold on