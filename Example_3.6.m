%% Minimal module basis/kernel representation from arbitrary basis (Example 3.6)
clear all;clc;
%------------------------Initialization------------------------------------
% G(z) is a given kernel representation/module basis
% G(z) = [z^2+1 (z+1)^2 z^3+z; z z+2 z^2];
%size of the matrices: G(z)-- 2*3 
p_G=2;q=3;N=3;L=4;l=1; %p_G=no of rows of G(z), q=no of cols of G(z) N=degree of G(z), L=length of the restricted behavior
% N+L = 7;
%------------------------------------------------------------------------

%create an array to initialize coefficients of G(z)
G0=[1 1 0; 0 2 0];G1=[0 2 1; 1 1 0];G2=[1 1 0; 0 0 1]; G3=[0 0 1;0 0 0]; %coefficient matrices of kernel rep
A0 = [G0 G1 G2 G3]; %first block of Toeplitz matrix

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

%-------------------Building left kernels (Algorithm 2)----------------------------
M_1=null(H_1')';%basis degree 0
s1=size(M_1);
bar_M_1=[zeros(s1(1,1),q) M_1;M_1 zeros(s1(1,1),q)];
M_2=null(H_2')';
s2=size(M_2);
basis_deg_1=[];
for j=1:s2(1,1)
    if(rank([bar_M_1;M_2(j,:)])== rank(bar_M_1))
        continue; 
    else basis_deg_1=[basis_deg_1;M_2(j,:)];
        bar_M_1=[bar_M_1;M_2(j,:)];
    end 
end
s3=size(bar_M_1);
bar_M_2=[zeros(s3(1,1),q) bar_M_1;bar_M_1 zeros(s3(1,1),q)];
M_3=null(H_3')';
s4=size(M_3);
basis_deg_2=[];
for j=1:s4(1,1)
    if(rank([bar_M_2;M_3(j,:)])== rank(bar_M_2))
        continue; 
    else basis_deg_2=[basis_deg_2;M_3(j,:)];
        bar_M_2=[bar_M_2;M_3(j,:)];
    end 
end
s5 = size(bar_M_2);
bar_M_3=[zeros(s5(1,1),q) bar_M_2;bar_M_2 zeros(s5(1,1),q)];
M_4=null(W')';
s6=size(M_4);
basis_deg_3=[];
for j=1:s6(1,1)
    if(rank([bar_M_3;M_4(j,:)])== rank(bar_M_3))
        continue; 
    else basis_deg_3=[basis_deg_3;M_4(j,:)];
        bar_M_3=[bar_M_3;M_4(j,:)];
    end 
end
% no degree 2 or degree 3 elements, only degree 0 and degree 1 annihilators 
%-------------------------------------------

%-----------minimal kernel representation--------------------
syms z;
r0 = M_1;
r1 = basis_deg_1(:,1:q)+z*basis_deg_1(:,q+1:2*q); 
min_ker = [r0; r1];
%-------------------------------------------------------------------------