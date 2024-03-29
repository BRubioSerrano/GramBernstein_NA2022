clear all
format rat
%Linear system Ax=b

%Experiments results of the
%Bidiagonal decomposition  of Gram matrix of 
%Bernstein  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 


%See experimental results in Mathematica: Gram_LinearSystem.nb

n=5;


A=zeros(n+1);
%Gram matrix of  Bernstein  basis  

for i=1:n+1
	for j=1:n+1
		A(i,j)=nchoosek(n,i-1)*nchoosek(n,j-1)*factorial(i+j-2)*factorial(2*n-i-j+2)/factorial(2*n+1); 
    end 
end

A
 
%Bidiagonal decomposition of Bernstein Gram matrix of 
%Geometric basis  

% BDA=BDAGram_matrix(n);
%  
% 
% %Linear system Ax=b
% 
% syms x
% 
% 
% v=1/(1+396*(x-0.5)^2);
% 
% x = -5:10/n:5
% 
% b=eval(v)
% 
% 
% %dlmwrite('b.csv',b,'precision','%.100f');
% 
% 
% % for i=1:n+1
% % fun1 = @(x) (1./(1+396.*(x-0.5).^2)).*nchoosek(n,i-1).*(1-x).^(n-i+1).*x.^(i-1);
% % b1(i) = integral(fun1,0,1);
% % end
% 
% for i=1:n+1
% fun1 = @(x) (0.01+x./(x.^2+1)).*nchoosek(n,i-1).*(1-x).^(n-i+1).*x.^(i-1);
% b1(i) = integral(fun1,0,1);
% end
% 
% 
% b=[30/10513,-63/12223,104/14095,-56/5891,89/7726,-87/6490,189/12473,-79/4715,25/1373,-167/8559,359/17371,-395/18222,...     
%     109/4834, -184/7901, 175/7321, -169/6925, 114/4597, -618/24625, 163/6441, -267/10496];
% % 
% % b=[1.40259225889142*10^-3, 2.460024678687843*10^-3, 3.500092610077788*10^-3, 4.517611751965842*10^-3, 5.507814767446761*10^-3, 6.466401859987952*10^-3,...
% % 7.389578218913938*10^-3, 8.274078391266732*10^-3, 9.117178065946993*10^-3, 9.916694118587380*10^-3, 1.067097405339338*10^-2, 1.137887618779081*10^-2,...
% % 1.203974205823599*10^-2, 1.265336258581963*10^-2,   1.321993953605909*10^-2,     1.374004374823445*10^-2,     1.421457150660234*10^-2,     1.464470028996585*10^-2,...
% % 1.503184497820909*10^-2, 1.537761542449372*10^-2,   1.568377612862638*10^-2,     1.595220857796437*10^-2,     1.618487666292518*10^-2,     1.638379542897204*10^-2,...
% % 1.655100329869643*10^-2,  1.668853778786975*10^-2, 1.679841464855291*10^-2,    1.688261030012201*10^-2,     1.694304735429596*10^-2,     1.698158300137060*10^-2]
% 
% 
% % size(b)
% 
% 
% % b=[2.853609003918695*10^(-3),    5.154214165485011*10^(-3),     7.378502675586492*10^(-3),     9.506023618475502*10^(-3),    1.151954186945111*10^(-2),...   
% %      1.340523268079367*10^(-2), 1.515272720173589*10^(-2), 1.675502629315182*10^(-2), 1.820830481499728*10^(-2),1.951163112405860*10^(-2),...
% %      2.066662704455418*10^(-2), 2.167709241279533*10^(-2), 2.254861586416712*10^(-2),2.328819025545957*10^(-2),2.390384740069584*10^(-2),...
% %      2.440432299104475*10^(-2), 2.479875894387648*10^(-2), 2.509644717921934*10^(-2), 2.530661607540821*10^(-2), 2.543825866509776*10^(-2)]
% % 
% 
%  %b= [1, 2, 3, 4, 1, 6, 7, 8, 9, 1, 6, 1, 2, 3, 3, 1,7,1,5,4];
% % 
% % 
% % 
% % % b=[6.359791964882841*10^(-4),8.486939356046927*10^(-4), 1.218518424570156*10^(-3),1.916953919719876*10^(-3),3.242501107004902*10^(-3),5.562459472462681*10^(-3),...
% % % 9.074806905203420*10^(-3), 1.345211561141169*10^(-2),1.766225780448371*10^(-2), 2.028759674181954*10^(-2), 2.028759674181954*10^(-2),1.766225780448371*10^(-2),...
% % % 1.345211561141169*10^(-2),9.074806905203417*10^(-3),5.562459472462680*10^(-3),3.242501107004901*10^(-3),1.916953919719876*10^(-3),1.218518424570156*10^(-3),...
% % % 8.486939356046925*10^(-4),6.359791964882843*10^(-4)];
% % 
% % 
% % 
% IB=TNInverseExpand(BDA);
% IM=inv(A);
% size(IB)
% size(b)
% SolB=IB*transpose(b)
% size(SolB)
% 
% 
% %  SolB=transpose(TNSolve(BDA,transpose(b)))
%    SolM=A\transpose(b)
% %   
%    dlmwrite('sistemaGramB.csv',SolB,'precision','%.100f');
%     dlmwrite('sistemaGramM.csv',SolM,'precision','%.100f');
% % % % % 
% % % %function TNSolve(B,b)
% % % %Solves a TN linear system Ax=b, where B=BD(A). (see TNSolve of Plamen Koev https://math.mit.edu/~plamen/software/TNTool.html)
% % % 
% % % %Using this bidiagonal decomposition, we can also obtaine the inverse, eigenvalues and singular values using the
% % % %functions presented in  https://math.mit.edu/~plamen/software/TNTool.html.
% % % 
% % % 
% % % 
% % % 
% % % 
