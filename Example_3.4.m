%% Example 3.4
clear all;clc;
syms z;
%a = 0.01; b = 0.05; eps = 0.001; %choice 1
a = 0.01; b = 0.05; eps = 0.0001; %choice 2
G = [(z+a)^2*(z+b)^2 0  -z z*(z+a)^2; -(z+a)^2*(z+b) z+b 0 -z^2] + eps*[-z^4 0 0 0;(-z^3+z^2) 0 0 0];
%------------Unimodular transformations on G(z)-----------
U1 = [(1+eps)/(1-eps) z;0 1];
G_1 = U1*G;
v1 = coeffs(G_1(1,1)); v2=coeffs(G_1(2,1));
U2 = [-v2(1,1)/v1(1,1) 1;0 1];
G_2 = U2*G_1; %1st row has degree 3, 2nd row has degree 3, checked by coeffs(G_2(i,j)).
%-------------------------------------------

%----------Data generation----------
G0 = [a^2*b^2 0 0 0;-a^2*b b 0 0]; G1 = [2*(a^2*b+a*b^2) 0 -1 a^2;-(a^2+2*a*b) 1 0 0];
G2=[a^2+4*a*b+b^2 0 0 2*a;-(2*a+b)+eps 0 0 -1]; G3 = [2*(a+b) 0 0 1;-(1+eps) 0 0 0]; 
G4= [1-eps 0 0 0;0 0 0 0];%coefficient matrices of kernel rep
A0 = [G0 G1 G2 G3 G4]; %first block of Toeplitz matrix

p_G=2;q=4;N=4;L=5; %p_G=no of rows of G(z), q=no of cols of G(z) N=degree of G(z), L=length of the restricted behavior

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

%----------------Check unimodularly transformed kernel rep on data----
P0 = subs(G_1, z, 0);
P1 = subs(diff(G_1), z, 0);
P2 = subs(diff(G_1,2)/2, z, 0);
P3 = subs(diff(G_1,3)/factorial(3), z, 0);
P4 = subs(diff(G_1,4)/factorial(4), z, 0);
P = [P0 P1 P2 P3];
P*W(1:(L-1)*q,:);%check if annihilates restricted behav 
%--------------------------------------------------------------

%--------------------Algorithm 2------------------------
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
M_3=null(H_3')';%basis degree 2, empty
M_4 = null(H_4',1e-10)';% basis degree 3
s4 = size(M_4);
%----------------------------------------------------------------------
