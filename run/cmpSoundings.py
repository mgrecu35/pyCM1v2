import matplotlib.pyplot as plt
import numpy as np
def read_s(fname):
    lines=open(fname,"r").readlines()
    vsL=[]
    for l in lines[1:]:
        vs=[float(v) for v in l.split()]
        vsL.append(vs[:3])
    return np.array(vsL),[float(v) for v in lines[0].split()]

s1,ss1=read_s("input_sounding.July23_99")
s2,ss2=read_s("input_sounding.June12_2011.CAPE_3000")

plt.plot(s2[:,2],np.array(s2[:,0]))
#plt.plot(s1[:,2],s1[:,0])
plt.ylim(0,1e4)

plt.figure()
q1int=np.interp(np.arange(180)*125+345,s1[:,0],s1[:,2])
q2int=np.interp(np.arange(180)*125+345,s2[:,0],s2[:,2])
t1int=np.interp(np.arange(180)*125+345,s1[:,0],s1[:,1])
t2int=np.interp(np.arange(180)*125+345,s2[:,0],s2[:,1])
plt.plot(q2int-q1int,np.arange(180)*125+500)

f=open("input_sounding.hyb","w")

f.write("%11.2f  %8.2f   %6.2f\n"%(ss2[0],ss2[1],ss2[2]))
for i in range(180):
    f.write("%11.2f  %8.2f   %6.2f   %6.2f   %6.2f\n"%\
            (345+i*125,0.5*(t1int[i]+t2int[i]),0.5*(q1int[i]+q2int[i]),0,0))

f.close()
 
