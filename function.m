function [ROI] = Emma_model2(pa,tt,func_hydrogel)
syms k1 k2 k3 k4 k5 k6 s kin C_h t
a = k5 + k1;
matrixA = [a+s, 0, 0, 0, 0, 0, 0;  k1, -(k2+s), 0, 0, 0, 0, 0;  0, k2, -(k3+s), 0, 0, 0, 0; 0, 0, k3, -(k4+s), 0, 0, 0;    0, 0, 0, k4, -s, 0, 0;    k5, 0, 0, 0, 0, -(k6+s), 0;    0, 0, 0, 0, 0, k6, -s];
matrixB = [kin*C_h; 0; 0; 0; 0; 0; 0];

% Cramer's Rule
matrix_v =  [matrixA(:,1:4), matrixB, matrixA(:,6:7)]; 
Cv_s = det(matrix_v)/det(matrixA);
%  Inverse Laplace Transform, from symbolic s to symbolic t
Cv_s = Cv_s/C_h;       Cv = ilaplace(Cv_s, s, t); 

% Substitute rate constants in
Cv_t1 = subs(Cv,[k1,k2,k3,k4,k5,k6,kin],pa);
C_vitreous = matlabFunction(Cv_t1);
t = tt;     % now set t as tt, the time vector input
ConcVit = C_vitreous(t); %  convert back from syms

ROI = conv(ConcVit,func_hydrogel); %  the solution is a convolution of the 
                                   %  "Drug input funtion" and ConcVit
