function[] = RK4_TEST_Plot(n, tmax, ic1, ic2)
    
    horiz = linspace(0, tmax, n);
    verti = RK4_TEST(n, tmax, ic1, ic2);
    
    plot(horiz, verti);

end