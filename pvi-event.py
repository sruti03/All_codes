import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math


magbx=np.load('nfmms1bx.30.npy')
magby=np.load('nfmms1by.30.npy')
magbz=np.load('nfmms1bz.30.npy')



time=np.load('Time.24.02.2019.mms1.npy')



l = time.size

#implementing the pvi technique
delBx = []
delBy = []
delBz = []

for i in range(l-32):
    delx = (magbx1[32+i]-magbx1[i])
    delBx = np.append(delBx,delx)
    dely = (magby1[32+i]-magby1[i])
    delBy = np.append(delBy,dely)
    delz = (magbz1[32+i]-magbz1[i])
    delBz = np.append(delBz,delz)


e=math.e

l1 = delBx.size

sumB = np.sqrt(np.sum((delBx)**2+(delBy)**2+(delBz)**2)/l1)           # PVI parameters
pvi = (np.sqrt((delBx)**2+(delBy)**2+(delBz)**2)/(sumB))
abspvi=abs(pvi)
p = pvi.size


pvi_mean = (np.sum(pvi))/l1                                          #mean PVI
std_pvi = np.sqrt(np.sum(((pvi)-(pvi_mean))**2)/l1)                  #standard deviation of PVI
kur_pvi = (np.sum(((pvi)-(pvi_mean))**4))/(((std_pvi)**4)*l1)        #kurtosis of PVI



#############################################################################

#eliminating events above threshold; checking the kurtosis
'''
c = (4.8*std_pvi)+pvi_mean
pvievent =[]

for i in range(l1):
    ele = pvi[i]
    if ele>=c:
        pvievent=np.append(pvievent,i)
index = pvievent.astype(int)
print(pvievent.size,l1)
#print(pvievent.size)
new_pvi = np.delete(pvi,index)
#print(new_pvi.size)
m = new_pvi.size

npvi_mean = (np.sum(new_pvi))/m
nstd_pvi = np.sqrt(np.sum(((new_pvi)-(npvi_mean))**2)/m)
nkur_pvi = (np.sum(((new_pvi)-(npvi_mean))**4))/(((nstd_pvi)**4)*m)

print(npvi_mean,nstd_pvi,nkur_pvi)
'''
#identifying events
pz=[]
qz=[]
srt_peakz=[]
end_peakz=[]

cz_max = (4*std_pvi)+pvi_mean
peaksz=[]
numz=0
for i in range(l1):
    if pvi[i]>=cz_max:
        numz=numz+1
    else:
        if numz>0:
            endz = i-1
            startz = endz-numz+1
            pz.append(startz)
            qz.append(endz)
        numz=0
s = len(pz)
m=[]
n=[]
for i in range(s):
    pks=np.arange(pz[i],qz[i]+1,1)
    u=pks.size
    v=u//2
    mks = pvi[pks]
    mk = np.argmax(mks)
    ind = pks[mk]
    m.append(ind)
    n=np.append(n,v)


for i in m:
    ct = i
    cz_min = (pvi[i])/e
    f=np.where(m==i)
    fi = int((n[f])/2)
    ks=i
    while ((pvi[ks-1]>cz_min) and ks>0):
        if ((pvi[ks-1]>pvi[ks]) and np.all(pvi[ks-fi:ks-1]>pvi[ks]) and (pvi[ks]<0.7*cz_max)):
            break
        else:
            ks=ks-1
    srt_peakz.append(ks)

for i in m:
    ct = i
    cz_min = (pvi[i])/e
    f=np.where(m==i)
    fi = int((n[f])/2)
    kf=i
    while ((pvi[kf+1]>cz_min) and kf<l1):
        if ((pvi[kf+1]>pvi[kf]) and np.all(pvi[kf+2:kf+fi+1]>pvi[kf]) and (pvi[kf]<0.7*cz_max)):
            break
        else:
            kf=kf+1
    end_peakz.append(kf)
    

#time duration of events
srt_time = time[srt_peakz]
end_time = time[end_peakz]
esize=len(srt_time)


time_durationz = time[end_peakz]-time[srt_peakz]

#plots
fig, ax = plt.subplots(1, 1,
                        figsize =(4, 3),
                        tight_layout = True)

plt.hist(time_durationz, bins = 20)
plt.xlabel('Duration')
plt.ylabel('Count per bin')
plt.title('MMS1 (30 rad/s) -- 1 sec')
plt.show()

