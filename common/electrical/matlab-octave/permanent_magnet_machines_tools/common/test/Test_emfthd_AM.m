
x = linspace(0,2,100);
y = 1*sin(x * pi);

test_slm = slmengine(x, y, 'Endcon', 'periodic');

thd = emfthd_AM(test_slm)

%%
x = linspace(0,2,100);
y = 1*sin(x * pi) + 0.1*sin(2*x * pi);

test_slm = slmengine(x, y, 'Endcon', 'periodic', 'knots', 10);

thd = emfthd_AM(test_slm)