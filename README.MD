# [iTrace](https://github.com/i5wear/iTrace)

高精度的要塞计算器，由 i5wear 制作。
An accurate stronghold calculator made by i5wear.

其研发基于要塞生成的所有细节及 [Cubiomes](https://github.com/Cubitect/Cubiomes)。
Its development is based on all details of stronghold generation and [Cubiomes](https://github.com/Cubitect/Cubiomes).

#### 使用帮助 User Guide

首次使用本程序时，请先设定游戏版本，并设定误差或开始误差校准。
If it's the first time to use, set the game version, then set error or start error calibration.

对准末影之眼按下 F3 + C，直接粘贴至本程序，即可快捷输入数据。
Aim at an ender eye and press F3 + C, then directly paste into the program for quick input.

若上述方式无法使用，则按下 F3，手动输入坐标和角度至本程序。
If the method is unavailable, press F3 and manually input positions and yaw into the program.

正常创建一个世界，获取种子并粘贴至本程序，即可开始误差校准。
Create a world as normal, get the seed and paste it into the program to begin calibration.

终止本程序不会丢失任何数据，所有历史输入均被自动保存与读取。
Exit won't cause any data loss, as all historical inputs are automatically saved and loaded.

#### 命令列表 Command List

添加末影之眼：ADD [X坐标] [Z坐标] [偏航角]
Add an ender eye: ADD [PosX] [PosZ] [Yaw]

设定游戏版本：VER [1.N]
Set game version: VER [1.N]

设定误差数值：ERR [平均数] [标准差]
Set error value: ERR [Mean] [SD]

开始误差校准：CAL [种子]
Begin calibration: CAL [Seed]

终止误差校准：CAL 0
End calibration: CAL 0

查看所有数据：CHECK
Check all data: CHECK

清除所有数据：CLEAR
Clear all data: CLEAR

终止程序运行：空输入
Exit program: Empty input
