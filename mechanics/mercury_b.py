

from pylab import*
import math

e = 0.2056
thet = 58.98/300000.0
theta = math.pow(thet,2)

time = 0.0
end_time = 100000.0
delta_t = 0.01
N = int(end_time/delta_t)
x = 1.0
y = 0.0
initial_v_x = 0.0
initial_v_y = 1.0

a_x = -x/(math.pow(x*x + y*y,1.5)*(1+e)) - (3*x*theta)/(math.pow(x*x + y*y,2.5)*(1+e))
a_y = -y/(math.pow(x*x + y*y,1.5)*(1+e)) - (3*y*theta)/(math.pow(x*x + y*y,2.5)*(1+e))
v_x = initial_v_x + a_x*delta_t/2.0
v_y = initial_v_y + a_y*delta_t/2.0

temp_x = x + v_x*delta_t/2.0
temp_y = y + v_y*delta_t/2.0

E =  (v_x*v_x + v_y*v_y)*(1+e)/2.0 - 1/math.sqrt(temp_x*temp_x + temp_y*temp_y)
L = math.fabs(initial_v_x*y - initial_v_y*x)
x_list = [x]
y_list = [y]
time_list = [time] 
energy_list = [E]
ang_list = [L]
for i in range(0,N):
    
    x += v_x*delta_t
    y += v_y*delta_t
    a_x = -x/(math.pow(x*x + y*y,1.5)*(1+e)) - (3*x*theta)/(math.pow(x*x + y*y,2.5)*(1+e))
    a_y = -y/(math.pow(x*x + y*y,1.5)*(1+e)) - (3*y*theta)/(math.pow(x*x + y*y,2.5)*(1+e))
    v_x += a_x*delta_t
    v_y += a_y*delta_t
    #temp_x = x + v_x*delta_t/2.0
    #temp_y = y + v_y*delta_t/2.0
    #E =  (v_x*v_x + v_y*v_y)*(1+e)/2.0 - 1/math.sqrt(temp_x*temp_x + temp_y*temp_y)
    #L = math.fabs(v_x*y - v_y*x)
    #time += delta_t
    x_list.append(x)
    y_list.append(y)
    #energy_list.append(E)
    #ang_list.append(L)
    #time_list.append(time)
    
a_list = x_list[0:10001]
b_list = y_list[0:10001]

c_list = x_list[4990001:5000001]
d_list = y_list[4990001:5000001]

e_list = x_list[9990001:10000001]
f_list = y_list[9990001:10000001]

sun_x = 0.0
sun_y = 0.0
plot(sun_x,sun_y,marker = 'o')
legend()
show()

xlabel('x')
ylabel('y')
title('Mercury\'s Orbbit')           
#plot(x_list,y_list,marker = '.',label = 'O')          
plot(a_list,b_list,marker = '.',label = 'O')
plot(c_list,d_list,marker = '.',label = 'O')
plot(e_list,f_list,marker = '.',label = 'O')
#plot(time_list,energy_list,marker = '.', label = 'E')
#plot(time_list,ang_list,marker = '.', label = 'L')
legend()
show()