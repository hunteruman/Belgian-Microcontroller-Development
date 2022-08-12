These results were used for the final correlation

This test includes:
Motor and cart for Ti = 0.05s and 0.1s
Inclination test (for Ti=0.05s): incl2 files

The most suspicious line of the code is:
volt_closed = feedback(sys_PI, sys_G);
So the simulated control signal differs from the measured one