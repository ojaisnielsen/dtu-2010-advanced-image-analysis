clear all
close all
clc
load pigsheart
colormap('gray')
imagesc(pigsheart);
hold on
% snake = FreeHandSnake();
snake = CircleSnake(size(pigsheart)/2+[-10,5], 500, min(size(pigsheart)/5));
DrawSnake(snake, 2);

%%
endsnake = evolveSnake(pigsheart, snake, 1, 3, 1, 0.2, 0.0001);
DrawSnake(endsnake, 3);
