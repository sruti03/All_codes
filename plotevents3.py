import numpy as np
import matplotlib.pyplot as plt

#Saving and plotting final events

time1=np.load('Time.24.02.2019.mms1.npy')
time2=np.load('Time.24.02.2019.mms2.npy')
time3=np.load('Time.24.02.2019.mms3.npy')
time4=np.load('Time.24.02.2019.mms4.npy')

deri1=np.load("nfderivativebx.15.mms1.npy")
deri2=np.load("nfderivativebx.15.mms2.npy")
deri3=np.load("nfderivativebx.15.mms3.npy")
deri4=np.load("nfderivativebx.15.mms4.npy")

l1=time1.size
l2=time2.size
l3=time3.size
l4=time4.size

deri_mean1 = (np.sum(deri1))/l1
std1 = np.sqrt(np.sum(((deri1)-(deri_mean1))**2)/l1)

deri_mean2 = (np.sum(deri2))/l2
std2 = np.sqrt(np.sum(((deri2)-(deri_mean2))**2)/l2)

deri_mean3 = (np.sum(deri3))/l3
std3 = np.sqrt(np.sum(((deri3)-(deri_mean3))**2)/l3)

deri_mean4 = (np.sum(deri4))/l4
std4 = np.sqrt(np.sum(((deri4)-(deri_mean4))**2)/l4)


sk1=np.load('peak1.npy')
sk2=np.load('peak2.npy')
sk3=np.load('peak3.npy')
sk4=np.load('peak4.npy')

ci1=np.load('2cI1.npy')
ci2=np.load('2cI2.npy')
ci3=np.load('2cI3.npy')
ci4=np.load('2cI4.npy')

fk1=np.load('endpeak1.npy')
fk2=np.load('endpeak2.npy')
fk3=np.load('endpeak3.npy')
fk4=np.load('endpeak4.npy')

hp1=np.load('hpeak1.npy')
hp2=np.load('hpeak2.npy')
hp3=np.load('hpeak3.npy')
hp4=np.load('hpeak4.npy')


pk1=np.delete(sk1,ci1)
pk2=np.delete(sk2,ci2)
pk3=np.delete(sk3,ci3)
pk4=np.delete(sk4,ci4)


ek1=np.delete(fk1,ci1)
ek2=np.delete(fk2,ci2)
ek3=np.delete(fk3,ci3)
ek4=np.delete(fk4,ci4)

hk1=np.delete(hp1,ci1)
hk2=np.delete(hp2,ci2)
hk3=np.delete(hp3,ci3)
hk4=np.delete(hp4,ci4)


srt_time1=time1[pk1]
srt_time2=time2[pk2]
srt_time3=time3[pk3]
srt_time4=time4[pk4]

end_time1=time1[ek1]
end_time2=time2[ek2]
end_time3=time3[ek3]
end_time4=time4[ek4]


der1=deri1[pk1]
der2=deri2[pk2]
der3=deri3[pk3]
der4=deri4[pk4]


esize1=srt_time1.size
esize2=srt_time2.size
esize3=srt_time3.size
esize4=srt_time4.size

# create 4 subplots
fig, ax = plt.subplots(nrows=2, ncols=2)

ax[0,0].plot(time1, deri1)
ax[0,1].plot(time2, deri2)
ax[1,0].plot(time3, deri3)
ax[1,1].plot(time4, deri4)

for i in range(esize1):
    w=pk1[i]
    if (deri1[w]>0):
        ax[0,0].axvline(x=srt_time1[i], ymin = 0.5, ymax = 1,color='orange')
        ax[0,0].axvline(x=hk1[i], ymin = 0.5, ymax = 1,color='green')
        ax[0,0].axvline(x=end_time1[i], ymin = 0.5, ymax = 1,color='orange')
    else:
        ax[0,0].axvline(x=srt_time1[i], ymin = 0, ymax = 0.5,color='orange')
        ax[0,0].axvline(x=hk1[i], ymin = 0, ymax = 0.5,color='green')
        ax[0,0].axvline(x=end_time1[i], ymin = 0, ymax = 0.5,color='orange')
ax[0,0].axhline(y=4*std1,color='red',ls='dashed',lw=0.5)
ax[0,0].axhline(y=-4*std1,color='red',ls='dashed',lw=0.5)
ax[0,0].axhline(y=0,color='black',ls='dashed',lw=0.5)

for i in range(esize2):
    w=pk2[i]
    if (deri2[w]>0):
        ax[0,1].axvline(x=srt_time2[i], ymin = 0.5, ymax = 1,color='orange')
        ax[0,1].axvline(x=hk2[i], ymin = 0.5, ymax = 1,color='green')
        ax[0,1].axvline(x=end_time2[i], ymin = 0.5, ymax = 1,color='orange')
    else:
        ax[0,1].axvline(x=srt_time2[i], ymin = 0, ymax = 0.5,color='orange')
        ax[0,1].axvline(x=hk2[i], ymin = 0, ymax = 0.5,color='green')
        ax[0,1].axvline(x=end_time2[i], ymin = 0, ymax = 0.5,color='orange')
ax[0,1].axhline(y=4*std2,color='red',ls='dashed',lw=0.5)
ax[0,1].axhline(y=-4*std2,color='red',ls='dashed',lw=0.5)
ax[0,1].axhline(y=0,color='black',ls='dashed',lw=0.5)

for i in range(esize3):
    w=pk3[i]
    if (deri3[w]>0):
        ax[1,0].axvline(x=srt_time3[i], ymin = 0.5, ymax = 1,color='orange')
        ax[1,0].axvline(x=hk3[i], ymin = 0.5, ymax = 1,color='green')
        ax[1,0].axvline(x=end_time3[i], ymin = 0.5, ymax = 1,color='orange')
    else:
        ax[1,0].axvline(x=srt_time3[i], ymin = 0, ymax = 0.5,color='orange')
        ax[1,0].axvline(x=hk3[i], ymin = 0, ymax = 0.5,color='green')
        ax[1,0].axvline(x=end_time3[i], ymin = 0, ymax = 0.5,color='orange')
ax[1,0].axhline(y=4*std3,color='red',ls='dashed',lw=0.5)
ax[1,0].axhline(y=-4*std3,color='red',ls='dashed',lw=0.5)
ax[1,0].axhline(y=0,color='black',ls='dashed',lw=0.5)

for i in range(esize4):
    w=pk4[i]
    if (deri4[w]>0):
        ax[1,1].axvline(x=srt_time4[i], ymin = 0.5, ymax = 1,color='orange')
        ax[1,1].axvline(x=hk4[i], ymin = 0.5, ymax = 1,color='green')
        ax[1,1].axvline(x=end_time4[i], ymin = 0.5, ymax = 1,color='orange')
    else:
        ax[1,1].axvline(x=srt_time4[i], ymin = 0, ymax = 0.5,color='orange')
        ax[1,1].axvline(x=hk4[i], ymin = 0, ymax = 0.5,color='green')
        ax[1,1].axvline(x=end_time4[i], ymin = 0, ymax = 0.5,color='orange')
ax[1,1].axhline(y=4*std4,color='red',ls='dashed',lw=0.5)
ax[1,1].axhline(y=-4*std4,color='red',ls='dashed',lw=0.5)
ax[1,1].axhline(y=0,color='black',ls='dashed',lw=0.5)


ax[0,0].set_title('Events:$B_x-MMS1$',fontsize=15)
ax[0,1].set_title('Events:$B_x-MMS2$',fontsize=15)
ax[1,0].set_title('Events:$B_x-MMS3$',fontsize=15)
ax[1,1].set_title('Events:$B_x-MMS4$',fontsize=15)

ax[1,0].set_xlabel('Time(secs)',fontsize=12)
ax[1,1].set_xlabel('Time(secs)',fontsize=12)

ax[0,0].set_ylabel(r'$\dfrac{dB_x}{dt}$',fontsize=12)
ax[1,0].set_ylabel(r'$\dfrac{dB_x}{dt}$',fontsize=12)

ax[0,0].set_xlim(7115,7125)
ax[0,0].set_ylim(-1.0,1.0)

ax[0,1].set_xlim(7115,7125)
ax[0,1].set_ylim(-1.0,1.0)

ax[1,0].set_xlim(7115,7125)
ax[1,0].set_ylim(-1.0,1.0)

ax[1,1].set_xlim(7115,7125)
ax[1,1].set_ylim(-1.0,1.0)

ax[0,0].text(7116, 0.8, '(a)', fontsize=10)
ax[0,1].text(7116, 0.8, '(b)', fontsize=10)
ax[1,0].text(7116, 0.8, '(c)', fontsize=10)
ax[1,1].text(7116, 0.8, '(d)', fontsize=10)


fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
fig.savefig('allmmsasi.eps',format='eps',dpi=300)
plt.show()

# saving final events
np.save('srt1',pk1)
np.save('srt2',pk2)
np.save('srt3',pk3)
np.save('srt4',pk4)
np.save('end1',ek1)
np.save('end2',ek2)
np.save('end3',ek3)
np.save('end4',ek4)


