function  test()
a1 = 1;
a2 = 1;
a3 = 0.1;
x = -a2:0.01:a2;
N1 = x.*(x - a2-a3)/((a1+a2)*(a1+2*a2+a3));
N2 = -(x + a1 + a2).*(x -a2 -a3)/((a1+a2)*(a2+a3));
N3 = x.*(x +a1 + a2)/((a2+a3)*(a1+2*a2+a3));
d1 = 5;
d2 = 3;
d3 = 3.5;
d = N1*d1+N2*d2+N3*d3;
plot(x,d);
end
