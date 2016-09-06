function [ROI] = Emma_model2(pa,tt,func_hydrogel,C_h)
syms k1 k2 k3 k4 k5 k6 s kin t
a = pa(5)+pa(1);

matrixA = [a+s, 0, 0, 0, 0, 0, 0;  pa(1), -(pa(2)+s), 0, 0, 0, 0, 0;  0, pa(2), -(pa(3)+s), 0, 0, 0, 0; 0, 0, pa(3), -(pa(4)+s), 0, 0, 0;    0, 0, 0, pa(4), -s, 0, 0;    pa(5), 0, 0, 0, 0, -(pa(6)+s), 0;    0, 0, 0, 0, 0, pa(6), -s]
matrixB = [C_h; 0; 0; 0; 0; 0; 0];

% Cramer's Rule
matrix_v =  [matrixA(:,1:4), matrixB, matrixA(:,6:7)]; 
Cv_s = det(matrix_v)/det(matrixA);
%  Inverse Laplace Transform, from symbolic s to symbolic t
Cv_s = Cv_s/C_h;       Cv = ilaplace(Cv_s, s, t) 

% Substitute rate constants in
C_vitreous = matlabFunction(Cv);
ConcVit = C_vitreous(tt) %  convert back from syms
func_hydrogel
ROI = transpose(func_hydrogel)*ConcVit; %  the solution is a convolution of the 
                                   %  "Drug input funtion" and ConcVit
