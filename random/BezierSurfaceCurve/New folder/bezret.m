function r = bezret(p)
n=length(p)-1;
for    i=0:1:n
sigma(i+1)=factorial(n)/(factorial(i)*factorial(n-i));  % for calculating (x!/(y!(x-y)!)) values 
end
l=[];
UB=[];
for u=0:0.002:1
for d=1:n+1
UB(d)=sigma(d)*((1-u)^(n+1-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end
P=l*p;
r=P;
%plot3(P(:,1),P(:,2),P(:,3),'b')
%hold on
%plot3(p(:,1),p(:,2),p(:,3),'b')