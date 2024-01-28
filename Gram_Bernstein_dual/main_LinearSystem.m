clear all

format longE
%Linear system Ax=b

%Experiments results of the
%Bidiagonal decomposition  of Gram matrix of dual  
%Bernstein  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 

%See experimental results in Mathematica: Gram_LinearSystem.nb

n=24
t=1;
l=2;

m=n-t-l+1


A=zeros(m+1);


%Gram matrix of  Bernstein  basis  

for i=1:m+1
	for j=1:m+1
		A(i,j)=nchoosek(n,t+i-1)*nchoosek(n,t+j-1)*factorial(2*t+i+j-2)*factorial(2*n-2*t-i-j+2)/factorial(2*n+1); 
    end 
end


 
% %Bidiagonal decomposition of Bersntein Dual Gram matrix of 
% %Geometric basis  
% 
BDA=BDAGram_Dual_matrix(n,m,t);

 
b=[17, -31, 77, -83, 27, -11, 96, -57, 70, -64, 29,-41,46,-16,74,-1,2,-6,7,-5,1,-2,6];


% 
% BDA
%Linear system Ax=b

% % 
% % % %  b=[17, -31, 77, -83, 27, -11, 96, -57, 70, -64, 29, -41,...
% % % %  46, -16, 74, -1, 2, -6, 7, -5, 1, -2, 6, -7, 5];
% % % 
 SolB=transpose(TNSolve(BDA,transpose(b)))
 SolM=A\transpose(b)

 
 
 
 dlmwrite('sistemaGramB.csv',SolB,'precision','%.45f');
  dlmwrite('sistemaGramM.csv',SolM,'precision','%.45f');
% % 
% % %function TNSolve(B,b)
% % %Solves a TN linear system Ax=b, where B=BD(A). (see TNSolve of Plamen Koev https://math.mit.edu/~plamen/software/TNTool.html)
% % 
% % %Using this bidiagonal decomposition, we can also obtaine the inverse, eigenvalues and singular values using the
% % %functions presented in  https://math.mit.edu/~plamen/software/TNTool.html.
% % 
% % 
% % 
% % 
% % 
