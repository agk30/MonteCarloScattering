from scipy.optimize import least_squares
import matplotlib.pyplot as plt

s_in = [0.1,3,1.9,3.2]
s_out = [1,2,3,4]

def residual(x, s_in, s_out):
    return s_out-(x*s_in)

answer = least_squares(residual,1,args=(s_in, s_out))

new_s_in = s_in*answer.x

plt.plot(s_in, 'r', s_out, 'b')

plt.plot(new_s_in, 'g', s_out, 'b')
plt.show()

print (answer.x)