# Polar-Codes
极化码是新近发现的理论上能达到香农限的唯一一种信道编码
0. 编码器采用蝶形图而非矩阵乘方法，相应模块代码在Polar_encoder.h/Polar_encoder.cpp;
1. 常用的译码算法有SC, SCL, CA-SCL，aCA-SCL译码算法，相应模块代码在Polar_decoder.h/Polar_decoder.cpp;
2. 常用的构造码字的方法有蒙特卡洛仿真、高斯近似等方法，相应模块代码在Polar_construction.h/Polar_construction.cpp;
3. 信道使用高斯信道，相应模块代码在Channel.h/Channel.cpp;
4. 调制方式为BPSK调制，相应模块代码在Modulation.h/Modulation.cpp;
5. 上层统计误码率、误比特率模块代码在ice_process.h/ice_process.cpp;
6. (1) main_TestConstruction.cpp 测试几种不同的构造方法;
   (2) main_TestDecodeSpeed.cpp 测试译码速度;
   (3) main_TestPerformaceCurve.CPP 统计译码曲线，可调用Matlab绘图接口显示误码率曲线。
