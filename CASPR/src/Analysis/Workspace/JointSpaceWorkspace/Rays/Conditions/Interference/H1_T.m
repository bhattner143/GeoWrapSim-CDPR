function out1 = H1_T(in1,in2,in3)
%H1_T
%    OUT1 = H1_T(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    14-Oct-2020 11:03:13

A1 = in2(1,:);
A2 = in2(2,:);
A3 = in2(3,:);
b11 = in1(1);
b12 = in1(4);
b21 = in1(2);
b22 = in1(5);
b31 = in1(3);
b32 = in1(6);
c1 = in3(1,:);
c2 = in3(2,:);
c3 = in3(3,:);
c4 = in3(4,:);
c5 = in3(5,:);
c6 = in3(6,:);
c7 = in3(7,:);
c8 = in3(8,:);
c9 = in3(9,:);
c10 = in3(10,:);
c11 = in3(11,:);
c12 = in3(12,:);
c13 = in3(13,:);
c14 = in3(14,:);
c15 = in3(15,:);
c16 = in3(16,:);
c17 = in3(17,:);
c18 = in3(18,:);
c19 = in3(19,:);
c20 = in3(20,:);
c21 = in3(21,:);
c22 = in3(22,:);
c23 = in3(23,:);
c24 = in3(24,:);
c25 = in3(25,:);
c26 = in3(26,:);
c27 = in3(27,:);
c28 = in3(28,:);
c29 = in3(29,:);
c30 = in3(30,:);
t2 = A1.^2;
t3 = A1.^3;
t4 = A2.^2;
t5 = A2.^3;
t6 = A3.^2;
t7 = A3.^3;
t8 = -b12;
t9 = -b22;
t10 = -b32;
t11 = A1+t8;
t12 = A2+t9;
t13 = A3+t10;
out1 = [b11.*c28+b21.*c29+b31.*c30+A1.*b11.*c22.*2.0+A2.*b11.*c25+A3.*b11.*c26+A2.*b21.*c23.*2.0+A1.*b21.*c25+A3.*b21.*c27+A1.*b31.*c26+A3.*b31.*c24.*2.0+A2.*b31.*c27+b11.*c1.*t3.*4.0+b11.*c6.*t5+b11.*c8.*t7+b11.*c13.*t2.*3.0+b21.*c2.*t5.*4.0+b21.*c4.*t3+b11.*c18.*t4+b11.*c20.*t6+b21.*c9.*t7+b21.*c14.*t4.*3.0+b21.*c16.*t2+b31.*c5.*t3+b31.*c3.*t7.*4.0+b31.*c7.*t5+b21.*c21.*t6+b31.*c17.*t2+b31.*c15.*t6.*3.0+b31.*c19.*t4+A1.*A2.*b11.*c16.*2.0+A1.*A3.*b11.*c17.*2.0+A1.*A2.*b21.*c18.*2.0+A2.*A3.*b21.*c19.*2.0+A1.*A3.*b31.*c20.*2.0+A2.*A3.*b31.*c21.*2.0+A2.*b11.*c4.*t2.*3.0+A3.*b11.*c5.*t2.*3.0+A1.*b11.*c10.*t4.*2.0+A1.*b11.*c11.*t6.*2.0+A1.*b21.*c6.*t4.*3.0+A2.*b21.*c10.*t2.*2.0+A3.*b21.*c7.*t4.*3.0+A2.*b21.*c12.*t6.*2.0+A1.*b31.*c8.*t6.*3.0+A3.*b31.*c11.*t2.*2.0+A2.*b31.*c9.*t6.*3.0+A3.*b31.*c12.*t4.*2.0,-c28.*t11-c29.*t12-c30.*t13-A1.*c22.*t11.*2.0-A2.*c23.*t12.*2.0-A1.*c25.*t12-A2.*c25.*t11-A1.*c26.*t13-A3.*c24.*t13.*2.0-A3.*c26.*t11-A2.*c27.*t13-A3.*c27.*t12-c1.*t3.*t11.*4.0-c2.*t5.*t12.*4.0-c4.*t3.*t12-c5.*t3.*t13-c6.*t5.*t11-c3.*t7.*t13.*4.0-c7.*t5.*t13-c8.*t7.*t11-c13.*t2.*t11.*3.0-c9.*t7.*t12-c14.*t4.*t12.*3.0-c16.*t2.*t12-c17.*t2.*t13-c18.*t4.*t11-c15.*t6.*t13.*3.0-c19.*t4.*t13-c20.*t6.*t11-c21.*t6.*t12-A1.*A2.*c16.*t11.*2.0-A1.*A3.*c17.*t11.*2.0-A1.*A2.*c18.*t12.*2.0-A2.*A3.*c19.*t12.*2.0-A1.*A3.*c20.*t13.*2.0-A2.*A3.*c21.*t13.*2.0-A2.*c4.*t2.*t11.*3.0-A3.*c5.*t2.*t11.*3.0-A1.*c6.*t4.*t12.*3.0-A1.*c10.*t4.*t11.*2.0-A2.*c10.*t2.*t12.*2.0-A3.*c7.*t4.*t12.*3.0-A1.*c8.*t6.*t13.*3.0-A1.*c11.*t6.*t11.*2.0-A3.*c11.*t2.*t13.*2.0-A2.*c9.*t6.*t13.*3.0-A2.*c12.*t6.*t12.*2.0-A3.*c12.*t4.*t13.*2.0];
