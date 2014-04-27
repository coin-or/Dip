close all
loc = [0, 2;
        0, 4;
        1, 2;
        1, 4;
        4, 1;
        4, 4;
        4, 5;
        5, 0;
        5, 2;
        5, 5];

plot(loc(:, 1), loc(:, 2), 'o')
hold on
for i = 1:size(loc, 1),
    text(loc(i, 1) + 0.1, loc(i, 2) + 0.1, int2str(i-1))
end;
axis([-0.5 5.5 -0.5 5.5])

figure

plot(loc(:, 1), loc(:, 2), 'o')
hold on
for i = 1:size(loc, 1),
    text(loc(i, 1) + 0.1, loc(i, 2) + 0.1, int2str(i-1))
end;
axis([-0.5 5.5 -0.5 5.5])

s1 = [0, 1, 3, 2];
s2 = [4, 7, 8];
s3 = [5, 6, 9];

plot([loc(s1 + 1, 1); loc(s1(1) + 1, 1)], [loc(s1 + 1, 2); loc(s1(1) + 1, 2)], 'c-.')
plot([loc(s2 + 1, 1); loc(s2(1) + 1, 1)], [loc(s2 + 1, 2); loc(s2(1) + 1, 2)], 'r--')
plot([loc(s3 + 1, 1); loc(s3(1) + 1, 1)], [loc(s3 + 1, 2); loc(s3(1) + 1, 2)], 'g:')

figure

plot(loc(:, 1), loc(:, 2), 'o')
hold on
for i = 1:size(loc, 1),
    text(loc(i, 1) + 0.1, loc(i, 2) + 0.1, int2str(i-1))
end;
axis([-0.5 5.5 -0.5 5.5])

s1 = [8, 9, 6, 5, 4, 7];

plot([loc(s1 + 1, 1); loc(s1(1) + 1, 1)], [loc(s1 + 1, 2); loc(s1(1) + 1, 2)], 'c-.')

figure

plot(loc(:, 1), loc(:, 2), 'o')
hold on
for i = 1:size(loc, 1),
    text(loc(i, 1) + 0.1, loc(i, 2) + 0.1, int2str(i-1))
end;
axis([-0.5 5.5 -0.5 5.5])

s1 = [5, 9, 6, 3, 1, 0, 2, 4, 7, 8];

plot([loc(s1 + 1, 1); loc(s1(1) + 1, 1)], [loc(s1 + 1, 2); loc(s1(1) + 1, 2)], 'b-')

