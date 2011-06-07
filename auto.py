import os

N       = 1001
epsilon = [-7.5+0.24j]
sigma   = [5e-9]
lambd   = [457.9e-7]
a       = [100e-9]
km      = [0]
kp      = [0]
km2     = [0]
kp2     = [0]
theta0  = [0., 20., 40., 60., 75.]
corr    = [1,2,3]
thetas  = [-90., 90.]

for i in theta0:
    File = open("par.in", "w")
    File.write("# First line\n")

    File.write(str(N) + "\n" \
               + str(i) + "\n" \
               + str(corr[0]) + "\n" \
               + str(a[0]) + "\n" \
               + str(km[0]) + "\n" \
               + str(kp[0]) + "\n" \
               + str(km2[0]) + "\n" \
               + str(kp2[0]) + "\n" \
               + str(epsilon[0].real) + "+i" + str(epsilon[0].imag) + "\n" \
               + str(sigma[0]) + "\n" \
               + str(thetas[0]) + "\n" \
               + str(thetas[1]) + "\n")
    File.close()
    os.system("nice -n 19 ./main")
    os.system("mv par.in  results/v" + str(i) + ".in")
    os.system("mv out.dat results/v" + str(i) + ".dat")
