%% Minimal module basis/kernel representation from arbitrary basis for MFD (Example 3.3)
clear all;clc;
%------------------------Initialization------------------------------------
% G(z) is a given kernel representation/MFD rep
% G(z) = [(z+1)^2(z+2)^2 0  -z z(z+1)^2; -(z+1)^2(z+2) z+2 0 -z^2];
%size of the matrices: G(z)-- 2*4 
p_G=2;q=4;N=4;L=5;l=3; %p_G=no of rows of G(z), q=no of cols of G(z) N=degree of G(z), L=length of the restricted behavior
% N+L = 9;
%------------------------------------------------------------------------
%create an array to initialize coefficients of G(z)
G0 = [4 0 0 0;-2 2 0 0]; G1 = [12 0 -1 1;-5 1 0 0]; G2=[13 0 0 2;-4 0 0 -1]; 
G3 = [6 0 0 1;-1 0 0 0]; G4= [1 0 0 0;0 0 0 0];%coefficient matrices of kernel rep
A0 = [G0 G1 G2 G3 G4]; %first block of Toeplitz matrix
%-----------------End of initialization----------------------------------------------------

%---------------------------Building toeplitz matrices iteratively--------

%Building toeplitz matrices Ai for frr rep G(z)
A = A0;
for j=1:N+L-1
    A = [A zeros(j*p_G,q); zeros(p_G,j*q) A0];
    %size(A)
end
%-----------------------End of toeplitz iterations--------------------------------------------------

%----------------Generating restricted behavior B_L (Algorithm 1)----------------------

Ker2 = null(A);% Kernels of non frr matrices at required stages
%W1=orth(Ker2);
W=orth(Ker2(1:L*q,:));
%--------------------------------------------------------------------------

%--------Building mosaic Hankel matrices --------------------------

%--------------Building mosaic H_1---------------
s1=size(W);
H_1 = [];
for i=1:s1(1,2)
    for k_1=0:L-1
        H_1 = [H_1 W(q*k_1+1:(k_1+1)*q,i)];
    end
end
%---------------------------------------------------

%--------------Building mosaic H_2-------------------
H_2=[];
for i=1:s1(1,2)
    for k_1=0:L-2
        H_2=[H_2 W(q*k_1+1:(k_1+2)*q,i)];
    end
end
%---------------------------------------------------

%--------------Building mosaic H_3-------------------
H_3=[];
for i=1:s1(1,2)
    for k_1=0:L-3
        H_3=[H_3 W(q*k_1+1:(k_1+3)*q,i)];
    end
end
%---------------------------------------------------

%--------------Building mosaic H_4-------------------
H_4=[];
for i=1:s1(1,2)
    for k_1=0:L-4
        H_4=[H_4 W(q*k_1+1:(k_1+4)*q,i)];
    end
end
%---------------------------------------------------

%-------------------Building left kernels (Algorithm 2)----------------------------
M_1 = null(H_1')';%basis degree 0, empty
M_2 = null(H_2')';%basis degree 1, empty
M_3=null(H_3')';%basis degree 2
s3 = size(M_3);
bar_M_3 = [zeros(s3(1,1),q) M_3;M_3 zeros(s3(1,1),q)];
M_4 = null(H_4')';
s4 = size(M_4);
basis_deg_3=[];
for j=1:s4(1,1)
    if(rank([bar_M_3;M_4(j,:)]) == rank(bar_M_3))
        continue; 
    else basis_deg_3 = [basis_deg_3;M_4(j,:)];
        bar_M_3 = [bar_M_3;M_4(j,:)];
    end 
end
s5 = size(bar_M_3);
bar_M_4 = [zeros(s5(1,1),q) bar_M_3;bar_M_3 zeros(s5(1,1),q)];
M_5 = null(W')';
s6 = size(M_5);
basis_deg_4 = [];
for j=1:s6(1,1)
    if(rank([bar_M_4;M_5(j,:)])== rank(bar_M_4))
        continue; 
    else basis_deg_4=[basis_deg_4;M_5(j,:)]; %no degree 4 annihilator
        %bar_M_4=[bar_M_3;M_4(j,:)];
    end 
end

%-----------minimal kernel representation--------------------
syms z;
r2 = M_3(:,1:q) + z*M_3(:,q+1:2*q) + z^2*M_3(:,2*q+1:3*q);
r3 = basis_deg_3(:,1:q) + z*basis_deg_3(:,q+1:2*q) + z^2*basis_deg_3(:,2*q+1:3*q) + z^3*basis_deg_3(:,3*q+1:4*q);
min_ker = [r2; r3];
s7 = size(basis_deg_3);
l1 = s3(1,2)/q-1; l2 = s7(1,2)/q-1;n=l1+l2; %l1,l2 are minimal lags, n = McMillan degree
%----------------------------------------------------------------

%--------------Unimpdular transformation method-----------------
G_z = [(z+1)^2*(z+2)^2 0  -z z*(z+1)^2; -(z+1)^2*(z+2) z+2 0 -z^2];
U1 = [1 z;0 1]; U2=[1 2; 0 1];
G_1z = U2*U1*G_z; %minimal kernel rep obtained via unimodular transformations
%-------------------------------------------------------------------------