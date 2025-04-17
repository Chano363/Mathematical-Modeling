clc;
clear;
close all;

S1 = readtable('data.xlsx','Sheet','问题2-曲线1');
t3_calculation(S1);
S2 = readtable('data.xlsx','Sheet','问题2-曲线2');
t3_calculation(S2);
S3 = readtable('data.xlsx','Sheet','问题2-曲线3');
t3_calculation(S3);