%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test of the function  trid_per_c2D              %
%  solving simultaneously m                        %
%   tridiagonal, periodic systems of size n        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
format long e;

n=input('Dimension of the matrix of the system n='); 
m=input('Number of sytems to solve             m=');

% Set the vectors defining the diagonals of matrices
for j=1:m
        % random values
    v1=rand(1,n); 
    v3=rand(1,n);
    v2=-(v1+v3);
        % coeff. of the matrix( 3 diagonals= 3 vectors)
    aa(j,:)=v1;
    ab(j,:)=v2;
    ac(j,:)=v3;
         % build the matrix of the system
    A=diag(v1(2:n),-1)+diag(v2,0)+diag(v3(1:n-1),1);
    A(1,n)=v1(1);
    A(n,1)=v3(n);
         % set the matrix diagonal dominant => invertible
    ab(j,:)=ab(j,:)+norm(A,inf)*ones(1,n);
    A=A+diag(norm(A,inf)*ones(1,n));
          % imposed "solution" (colum vectors)
    sol(:,j)=[j:j+n-1]';
          % corrsponding RHS term
    fi(:,j)=A*sol(:,j);
          % Matlab solution
    solmat(:,j)=A\fi(:,j);
end

% Solution using trid_per_2
soltrid=NSE_trid_per_c2D(aa,ab,ac,fi');

% Compare the solutions
for j=1:m
   fprintf('===== System  number  %i\n',j)
   fprintf('trid_per_c2D //   Matlab //   Exact\n')
    disp([soltrid(j,:)' solmat(:,j) sol(:,j)])
  fprintf('Press return to continue\n');pause
end