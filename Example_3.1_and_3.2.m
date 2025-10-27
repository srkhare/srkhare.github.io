%% Example 3.1 and 3.2
clear all;clc;
%------------------------Initialization------------------------------------
%G(z) = [z^2 0 -1 0; 0 z 0 -1; z^3 z -z -1];
%size of G(z)-- 3*4, lag = 2 
p_G=3;q=4;N=3;L=4; %p_G=no of rows of G(z), q=no of cols of G(z) N=degree of G(z), L=length of the restricted behavior
p_R=2;l=2;%p_R=no of rows of R(z), l=lag
%N+L = 7;
%------------------------------------------------------------------------

%create an array to initialize coefficients of G(z)
G0=[0 0 -1 0;0 0 0 -1;0 0 0 -1];G1=[0 0 0 0; 0 1 0 0 ;0 1 -1 0];G2=[1 0 0 0;0 0 0 0;0 0 0 0];G3=[0 0 0 0;0 0 0 0;1 0 0 0];%coefficient matrices of non frr kernel rep
% G_0..G_N are stored in G(1)...G(N+1)
A0=[G0,G1,G2,G3];
%-----------------End of initialization----------------------------------------------------

%---------------------------Building toeplitz matrices iteratively--------

%Building toeplitz matrices Ai for frr rep G(z)
A = A0;
for j=1:L-1
    A = [A zeros(j*p_G,q); zeros(p_G,j*q) A0];
    %size(A)
end
%-----------------------End of toeplitz iterations--------------------------------------------------

%----------------Generating restricted behavior B_L (Algorithm 1)----------------------

Ker2 = null(A);% Kernels of non frr matrices at required stages
%W1=orth(Ker2);
W = orth(Ker2(1:L*q,:)); % orthonormal basis of the restricted behavior
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

%-------------------Building left kernels----------------------------
M_2=null(H_2')';%basis degree 1
s2=size(M_2);
bar_M_2=[zeros(s2(1,1),q) M_2;M_2 zeros(s2(1,1),q)];
M_3=null(H_3')';
s3=size(M_3);
basis_deg_2=[];
for j=1:s3(1,1)
    if(rank([bar_M_2;M_3(j,:)])== rank(bar_M_2))
        continue; 
    else basis_deg_2=[basis_deg_2;M_3(j,:)];
        bar_M_2=[bar_M_2;M_3(j,:)];
    end 
end
bar_M_3=[zeros(s3(1,1),q) M_3;M_3 zeros(s3(1,1),q)];
%bar_M_3 is 6*16, W is 16*%, bar_M_3*W =0 verified.
%---------------------------------

%-----------minimal kernel representation--------------------
syms z;
r1 = M_2(:,1:q) + z*M_2(:,q+1:2*q);
r2 = basis_deg_2(:,1:q) + z*basis_deg_2(:,q+1:2*q) + z^2*basis_deg_2(:,2*q+1:3*q);
min_ker = [r1; r2];
s4 = size(basis_deg_2);
l1 = s2(1,2)/q-1; l2 = s4(1,2)/q-1;n=l1+l2; %l1,l2 are minimal lags, n = McMillan degree
%----------------------------------------------------------------