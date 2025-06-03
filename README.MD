# [iTrace](https://github.com/i5wear/iTrace)

高精度的要塞计算器，改进了 Ninjabrain Bot 的算法。

其研发基于要塞生成的所有细节，它们均来自 [Cubiomes](https://github.com/Cubitect/Cubiomes)。

An accurate stronghold calculator, which improved the algorithm of Ninjabrain Bot.

Its development is based on all details of stronghold generation that come from [Cubiomes](https://github.com/Cubitect/Cubiomes).

#### 使用帮助 User Guide

首次使用时，请先设定游戏版本，并设定误差或开始误差校准。

对准末影之眼按下 F3 + C，直接将其粘贴，即可快捷输入数据。

若上述方式无法使用，则应按下 F3，并手动输入坐标和角度。

正常创建一个世界，获取种子并将其粘贴，即可开始误差校准。

退出程序不会丢失任何数据，所有输入均被自动保存与读取。

At the first time to use, set the game version, then set error or start calibration.

Aim at an ender eye and press F3 + C, then paste it directly for quick data input.

If this method is unavailable, press F3 and input positions and angle manuallly.

Create a world as normal, get the seed and paste it to begin error calibration.

Exit program does not cause any data loss, as all inputs are auto saved and loaded.

#### 命令列表 Command List

添加末影之眼：ADD [X坐标] [Z坐标] [偏航角]

设定游戏版本：VER [1.N]

设定误差数值：ERR [平均数] [标准差]

开始误差校准：CAL [种子]

终止误差校准：CAL 0

查看所有数据：CHECK

清除所有数据：CLEAR

终止程序运行：空输入

Add an ender eye: ADD [PosX] [PosZ] [Yaw]

Set game version: VER [1.N]

Set error value: ERR [Mean] [SD]

Begin calibration: CAL [Seed]

End calibration: CAL 0

Check all data: CHECK

Clear all data: CLEAR

Exit program: Empty input
