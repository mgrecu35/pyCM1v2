from numpy import *

q2_MCS=[2.36484,-955.979,4.72795,-930.660,7.09086,-907.678,9.78446,-870.725,11.5981,-849.989,12.4811,-801.046,11.2782,-747.079,10.0215,-681.415,9.42340,-618.199,8.22311,-536.181,5.81737,-428.247,2.80878,-308.524,0.0163627,-226.240]

q1_MCS=[10.1450,-948.045,-1.48485,-861.175,-8.55755,-767.034,-14.5838,-662.326,-22.3911,-512.196,-26.2998,-428.377,-28.1096,-333.929,-22.5429,-253.058,-14.4976,-210.565,-6.79707,-178.597,-0.859038,-132.725]

ls=open("q1q2.sesa","r").readlines()
q1q2L=[]
for i in range(244):
    dsnd=[]
    for l1 in ls[i*42+1:i*42+42]:
        v3=[float(v) for v in l1.split()[:3]]
        dsnd.append(v3)
    q1q2L.append(dsnd)

