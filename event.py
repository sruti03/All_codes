import numpy as np
import matplotlib.pyplot as plt
import math

e=math.e
mag=np.load('nfsmooth.magneticbx.15.mms4.npy')              # smooth magnetic field
deri_bx = np.load('nfderivativebx.15.mms4.npy')             # derivative of magnetic field

time = np.load('Time.24.02.2019.mms4.npy')                  # time

l = time.size

#mean, standard deviation of magnetic field
deribx_mean = (np.sum(deri_bx))/l
abs_deribx = abs(deri_bx)
std_devbx = np.sqrt(np.sum(((deri_bx)-(deribx_mean))**2)/l)

##################
#events present above 4sigma threshold value

pz=[]
qz=[]
srt_peak=[]         #start of event index
end_peak=[]         #end of event index

c_max = 4*std_devbx

peaks=[]
num=0
for i in range(l):
    if abs_deribx[i]>=c_max:
        num=num+1
    else:
        if num>0:
            end = i-1
            start = end-num+1
            pz.append(start)
            qz.append(end)
        num=0
s = len(pz)
m=[]
n=[]
for i in range(s):
    pks=np.arange(pz[i],qz[i]+1,1)
    u=pks.size
    v=u//2
    mks = abs_deribx[pks]
    mk = np.argmax(mks)
    ind = pks[mk]
    m.append(ind)
    n=np.append(n,v)

for i in m:
    ct = i
    c_min = (abs_deribx[i])/e
    f=np.where(m==i)
    fi = int((n[f])/2)
    ks = i
    while ((abs_deribx[ks-1]>c_min) and ks>0):
        if ((abs_deribx[ks-1]>abs_deribx[ks]) and np.all(abs_deribx[ks-fi:ks-1]>abs_deribx[ks]) and (abs_deribx[ks]<0.7*c_max)):
            break
        else:
            ks=ks-1
    srt_peakz.append(ks)

for i in m:
    ct = i
    c_min = (abs_deribx[i])/e
    f=np.where(m==i)
    fi = int((n[f])/2)
    kf = i
    while ((abs_deribx[kf+1]>c_min) and kf<l):
        if ((abs_deribx[kf+1]>abs_deribx[kf]) and np.all(abs_deribx[kf+2:kf+fi+1]>abs_deribx[kf]) and (abs_deribx[kf]<0.7*c_max)):
            break
        else:
            kf=kf+1
    end_peakz.append(kf)

srt_time = time[srt_peakz]                 #start of event
end_time = time[end_peakz]                 #end of event
esize=len(srt_time)
highpeak=time[m]                           #peak of an event

np.save("peak4",srt_peakz)
np.save("endpeak4",end_peakz)
np.save("hpeak4",highpeak)