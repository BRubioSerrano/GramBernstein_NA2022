function BDA =BDAGram_Dual_matrix(n,m,t)

%Bidiagonal decomposition of the Gram matrix of 
%Bernstein basis 

BDA=zeros(m+1);


%Computation of the multipliers m_{i,j}

 for i=2:m+1
      for j=1:i
          BDA(i,j)=(2*t+i-1)*(n-t-i+2)*(2*n-2*t-i+3)/((t+i-1)*(2*n-2*t-i-j+3)*(2*n-2*t-i-j+4));
      end
 end  

%Computation of the pivots p_{i,i}
 for i=1:m+1
     BDA(i,i)=(nchoosek(n,t+i-1)*nchoosek(n,t+i-1)*factorial(2*t+i-1)*factorial(i-1)*factorial(2*n-2*t-2*i+2)*factorial(2*n-2*t-2*i+3))/(factorial(2*n-i+2)*factorial(2*n-2*t-i+2));
 end
 
 

%Computation of the multipliers tilde m_{i,i}
 
for i=1:m+1
    for j=i+1:m+1
          BDA(i,j)=(2*t+j-1)*(n-t-j+2)*(2*n-2*t-j+3)/((t+j-1)*(2*n-2*t-i-j+3)*(2*n-2*t-i-j+4)); 
     end
end


