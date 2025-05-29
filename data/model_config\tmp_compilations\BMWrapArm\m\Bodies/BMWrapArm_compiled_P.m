function out1 = BMWrapArm_compiled_P(in1,in2,in3,in4)
out1 = zeros(6,6);
out1(1:6, 1:6) = BMWrapArm_compiled_P_1_1(in1,in2,in3,in4);
