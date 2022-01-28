import numpy as np

ls=open("OK_sounding_July23_99").readlines()
ls=open("OK_sounding_June12_2011").readlines()

vs=[float(v) for v in ls[0].split()]
snd_dataL=[]
cmp_snd=[]
import metpy.calc as mpcalc
from metpy.units import units

for i,l in enumerate(ls[1:-2]):
    v=[float(v) for v in l.split()]
    if i>37 and i<53:
        print(l)
        f=(481-v[0])/(481-269)
        rh_int=(1-f)*71+f*30
        ws=mpcalc.mixing_ratio_from_relative_humidity(v[0]*units.hectopascal,\
                                                      v[2]*units.celsius,\
                                                      rh_int*units.percent)
        print(rh_int,ws)
        ws=ws.magnitude*1e3
    else:
        ws=v[5]
    snd_dataL.append([v[1],v[8],ws,0.0,0.0])
    cmp_snd.append([v[0],v[2],v[3],v[4],v[1]])

cmp_snd=np.array(cmp_snd)
from tephigram import Tephigram

tephigram = Tephigram()


tephigram.plot_sounding(P=cmp_snd[:,0],T=cmp_snd[:,1],T_dp=cmp_snd[:,2])
parcel_info = tephigram.plot_test_parcel(z=cmp_snd[:,-1], P=cmp_snd[:,0], \
                                         T=cmp_snd[:,1], RH=cmp_snd[:,3]/100.)
tephigram.plot_legend()
f=open("input_sounding","w")
#stop
f.write("%11.2f  %8.2f   %6.2f\n"%(vs[0],snd_dataL[0][1],snd_dataL[0][2]))
for i in range(len(snd_dataL)):
    if snd_dataL[i][1]<3000:
        f.write("%11.2f  %8.2f   %6.2f   %6.2f   %6.2f\n"%\
                (snd_dataL[i][0],snd_dataL[i][1],1.0*snd_dataL[i][2],0,0))
    else:
        if snd_dataL[i][1]<6000:
            f.write("%11.2f  %8.2f   %6.2f   %6.2f   %6.2f\n"%\
                    (snd_dataL[i][0],snd_dataL[i][1],1.0*snd_dataL[i][2],0,0))
        else:
            f.write("%11.2f  %8.2f   %6.2f   %6.2f   %6.2f\n"%\
                    (snd_dataL[i][0],snd_dataL[i][1],1.0*snd_dataL[i][2],0,0))

f.close()
       
