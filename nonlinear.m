%% BFSfinder
% 3.4 start
clear
clc
A1=[1 1 -1 0 0 0 0];
A2=[1 -1 0 1 0 0 0];
A3=[-1 1 0 0 1 0 0];
A4=[1 0 0 0 0 1 0];
A5=[0 1 0 0 0 0 1];
A=[A1;A2;A3;A4;A5];
B1=[A(:,1) A(:,2) A(:,3) A(:,4) A(:,5)] %% a basis
b=[1;1;1;2;2]
x1=B1\b  %x1=2 x2=2 s1=3 s2=1 s3=1 is the solition to 3.4
%% 3.4 continued
%remove a basis, let's say s2 and add d5
B2=B1;

B2(:,4)=A(:,7) % vi har x1 x2 s1 s5 s3
A
x2=B2\b % new bfs

B3=B1;
B3(:,5)=A(:,6) % remove s3, add s4
x3=B3\b 

% [B3 A(:,5) A(:,7)]*[x3;0;0] %% it works!
% %% tjafs
% B2*x2
% B3=B2;
% B3(:,3)=A(:,6) % vi har bytt till x1 x2 s4 s5 s3
% x3=B3\b % en ny bfs är funnen
% B3=B1;
% B3(:,2)=A(:,2);
% B3=[A(:,1) A(:,2) A(:,4) A(:,5) A(:,7)]; % not feasible
% B4=[A(:,3) A(:,4) A(:,5) A(:,6) A(:,7)];
% x3=B4\b;

%% 9.1
A1=[1;0];
A2=[0;1];
X1=[-3;1];
X2=[-2;1];
X3=[1;-2];
S1=[-1;0];
S2=[0;-1];
b=[3;1]
% 
B=[A1 A2]
N=[X1 X2 X3 S1 S2]
cB=[1;1]
cN=[0;0;0;0;0]

CtildeN=cN'-cB'*B^-1*N  % we se all ctilde_n>0--> optimal solution xb=
xb=B\b
%% 9.1
A1=[1;0];
A2=[0;1];
X1=[-1;-2];
X1neg=[1;2];
X3=[2;-2];
S1=[-1;0];
S2=[0;-1];
b=[0.5;1]
% 
B=[A1 A2]
N=[X1 X1neg X3 S1 S2]
cB=[1;1]
cN=[0;0;0;0;0]

CtildeN=cN'-cB'*B^-1*N  % we se ctilde_n(2)<0--> optimal solution xb not found, we now update!
% we know that X1neg should be swithced, but to what?
% check
a=B^-1*b % 0.5 1 %larger than 0
c=B\X1neg
% 
%% 

   % (B\b)%% /(B\N)(i)

a(1)/c(1) % 0.5, check the other
a(2)/c(2) %=0.5 so we may choose either a1 or a2 to leave the model,
% switch a1 to x1neg
% New guess
B=[X1neg A2]
N=[X1 A1 X3 S1 S2]
cB=[0;1];
cN=[0;1;0;0;0];
CtildeN=cN'-cB'*B^-1*N % we get S1 negative, it should enter the basis,
c=B\S1% =-1 2, whch means only base element 2 can leave the basis since the other is negatve
% so: a2 is replaced by s1, new base becomes x1neg,s1. now: since a1 and a2
% are not in the basis, they must be= 0,0, so we found a BFS! we are now
% done with phase 1. 

%% PHASE 2 9.4
% start with base x1neg, s1
B=[X1neg S1]
N=[X1 X3 S2]
cB=[1;0]
cN=[-1;1;0]; % this stuff relates to how it looks in the minimization function! 
CtildeN=cN'-cB'*B^-1*N %>or equal to 0 for every solution! so we are at optimal solution!
%%
clc
clear
A1=[1;0];
A2=[0;1];
X1=[-1;-0];
X2=[3;1];
S1=[-1;0];
S2=[0;1];
b=[3;2]
B=[A1 A2]
N=[X1 X2 S1 S2]
cB=[1;1];
cN=[0;0;0;0];
CtildeN=cN'-cB'*B^-1*N; % 1    -4     1    -1 --> x2 is incoming variable!
% calculate B^-1*b to check larger than zero!
B\b
(B\b)./(B\X2) % now take minimum value, this corresponds to outgoing base! 1,2 . 1 is smallest, A1 goes out!
%% New BFS!
B=[X2 A2]
N=[X1 A1 S1 S2]
cB=[0;1] % cB-[0;1]
cN=[0;1;0;0] % cN+[0;1;0;0] 
CtildeN=cN'-cB'*B^-1*N % -0.3333    1.3333   -0.3333   -1.0000 % smallest at s2! s2 icoming
B\b % above zero
%%
(B\b)./(B\S2) % infinite for first, so we take second! means we msut switch out a2! new
%% 4.5 PHASE 2
X1=[-1;0];
X2=[3;1];
S1=[-1;0];
S2=[0;1];
B=[X2,S2]
N=[X1 S1]
cB=[1;0];  %x2 is in real function, s2 is not 
cN=[-2;0]; % x1 is in real, s1 is not
CtildeN=cN'-cB'*B^-1*N %=-1.667 0.33. less than zero for j=1, x1 is incoming variable!
(B\b)./(B\X1) % = -3 3--> x2 will be outgoing variable!
%% continued
B=[X1,S2]
N=[X2 S1]
cB=[-2;0];  %x2 is in real function, s2 is not 
cN=[1;0];
CtildeN=cN'-cB'*B^-1*N % -5 2, x2 is incoming variable
(B\b) % larger than 0 for second entry, s2 bust be outgoing,
%% continued to optimality
% Now we have a new base! X2 and S2 switched!
B=[X1,X2]
N=[S2 S1]
cB=[-2;1]
cN=[0;0]
CtildeN=cN'-cB'*B^-1*N %% now ctilde 5 2 is larger than 0 for borth values--> optimal solution!
xoptimal=B\b
xoptimal'*[-2 1]'
%% 4.12
clc
clear
b=[2;1]
A1=[1;0];
A2=[0;1];
X1pos=[-1;-3];
X1neg=[1;3];
X2=[-2;-1];
X3=[-1;0];
S1=[0;-1]
B=[A1 A2]
N=[X1pos X1neg X2 X3 S1]
cB=[1;1];
cN=[0;0;0;0;0]
CtildeN=cN'-cB'*B^-1*N % 4    -4     3     1     1 --> x1neg is incoming
%B\b is positive, so we have a bounded direction, chose outgoing variable
(B\b)./(B\X1neg) % =2, 0.333 so a2 must be outgoing
%% continued new base
clc
B=[A1 X1neg]
N=[X1pos A2 X2 X3 S1]
B\b
cB=[1;0];% x1neg is 0 in obj function
cN=[0;1;0;0;0];% a2 is in second position
CtildeN=cN'-cB'*B^-1*N % ony negative for S1, it is incoming
(B\b)./(B\S1) %  negative for i = 2, when B\b is positive, so A1 is outgoing variable
%% 
B=[S1 X1neg];
N=[X1pos A2 X2 X3 A1];
B\b
cB=[0;0];% x1neg is 0 in obj function
cN=[0;1;0;0;1] % a2 is in second position
CtildeN=cN'-cB'*B^-1*N % only positive! we have optimal solution! --> bfs found
%% PHASE 2 4.12
clc
clear
b=[2;1]
X1pos=[-1;-3];
X1neg=[1;3];
X2=[-2;-1];
X3=[-1;0];
S1=[0;-1];
B=[S1 X1neg];
N=[X1pos X2 X3];
cB=[0;1];% x1neg is 1, s1 is 0
cN=[-1;-1;0] %  x1pos and x2 are negative in obj functin is in second position
CtildeN=cN'-cB'*B^-1*N % same in i= 1 and 2, but non-negative! optimal solution found! still chose x1pos as incoming
B\b % 5 2 >0, so it is unique! 
[0;1]'*(B\b) % optimal value=2
(B\b)./(B\X1pos) % not possible to get new point,
%% get new set!