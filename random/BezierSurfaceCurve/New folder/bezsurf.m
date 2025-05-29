P00=[0 0 0];
P01=[-5 18 5];
P02=[-7 32 7];
P03=[0 50 -3];
P10=[12 0 -5];
P11=[13 17 -3];
P12=[14 31 10];
P13=[12 55 5];
P20=[28 1 4];
P21=[29 15 3];
P22=[28 28 0];
P23=[29 45 -7];
P30=[40 0 1];
P31=[43 18 -1];
P32=[37 32 5];
P33=[40 50 0];
p0=[P00;P01;P02;P03];
p1=[P10;P11;P12;P13];
p2=[P20;P21;P22;P23];
p3=[P30;P31;P32;P33];
%write 4 curves
r0=bezret(p0);
r1=bezret(p1);
r2=bezret(p2);
r3=bezret(p3);
B=[];
n=length(r0);
for i=1:n
    B(1:4,(3*i-2):3*i)=[r0(i,:);r1(i,:);r2(i,:);r3(i,:)];
    bez3d(B(1:4,(3*i-2):3*i))
end


%{
R=bez3d(p);
R=bez3d(p);
R=bez3d(p);
%arrange lementwise
%run bezier on all and plot
%}