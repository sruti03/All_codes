import numpy as np
from itertools import chain

#Eliminating sporadic events by different criteria

time1=np.load('Time.24.02.2019.mms1.npy')
time2=np.load('Time.24.02.2019.mms2.npy')
time3=np.load('Time.24.02.2019.mms3.npy')
time4=np.load('Time.24.02.2019.mms4.npy')

deri1=np.load("nfderivativebx.15.mms1.npy")
deri2=np.load("nfderivativebx.15.mms2.npy")
deri3=np.load("nfderivativebx.15.mms3.npy")
deri4=np.load("nfderivativebx.15.mms4.npy")

pk1=np.load('peak1.npy')
pk2=np.load('peak2.npy')
pk3=np.load('peak3.npy')
pk4=np.load('peak4.npy')

smms1=time1[pk1]
smms2=time2[pk2]
smms3=time3[pk3]
smms4=time4[pk4]

der1=deri1[pk1]
der2=deri2[pk2]
der3=deri3[pk3]
der4=deri4[pk4]

l1=smms1.size
l2=smms2.size
l3=smms3.size
l4=smms4.size


# distance between different satellites
solar_speed = 330
m23 = 180/solar_speed
m14 = 30/solar_speed
m12 = 60/solar_speed
m43 = 90/solar_speed
m13 = 120/solar_speed
m24 = 90/solar_speed


#eliminating events w.r.t. mms1
c1 = []
for i in range(l1):
    ele= smms1[i]
    eer2 = (np.array(np.where(np.logical_and(smms2>=ele-m12, smms2<=ele+m12))))
    eer3 = (np.array(np.where(np.logical_and(smms3>=ele-m13, smms3<=ele+m13))))
    eer4 = (np.array(np.where(np.logical_and(smms4>=ele-m14, smms4<=ele+m14))))
    er2 = list(chain(*eer2))
    er3 = list(chain(*eer3))
    er4 = list(chain(*eer4))
    rc2=[]
    rc3=[]
    rc4=[]
    r2=0
    if eer2.size!=0:
        for j1 in range(eer2.size):
            j=er2[j1]
            if (der1[i]<0 and der2[j]<0):
                rc2.append(j)
                r2=len(rc2)
            elif (der1[i]>0 and der2[j]>0):
                rc2.append(j)
                r2=len(rc2)
            else:
                r2=0

    r3=0
    if eer3.size!=0:
        for j1 in range(eer3.size):
            j=er3[j1]
            if (der1[i]<0 and der3[j]<0):
                rc3.append(j)
                r3=len(rc3)
            elif (der1[i]>0 and der3[j]>0):
                rc3.append(j)
                r3=len(rc3)
            else:
                r3=0

    r4=0
    if eer4.size!=0:
        for j1 in range(eer4.size):
            j=er4[j1]
            if (der1[i]<0 and der4[j]<0):
                rc4.append(j)
                r4=len(rc4)
            elif (der1[i]>0 and der4[j]>0):
                rc4.append(j)
                r4=len(rc4)
            else:
                r4=0

    if (r2==0 and r3==0) or (r3==0 and r4==0) or (r2==0 and r4==0):
        c1.append(i)


##eliminating events w.r.t. mms2
c2 = []
for i in range(l2):
    ele= smms2[i]
    eer1 = (np.array(np.where(np.logical_and(smms1>=ele-m12, smms1<=ele+m12))))
    eer3 = (np.array(np.where(np.logical_and(smms3>=ele-m23, smms3<=ele+m23))))
    eer4 = (np.array(np.where(np.logical_and(smms4>=ele-m24, smms4<=ele+m24))))
    er1 = list(chain(*eer1))
    er3 = list(chain(*eer3))
    er4 = list(chain(*eer4))
    rc1=[]
    rc3=[]
    rc4=[]
    r1=0
    if eer1.size!=0:
        for j1 in range(eer1.size):
            j=er1[j1]
            if (der2[i]<0 and der1[j]<0):
                rc1.append(j)
                r1=len(rc1)
            elif (der2[i]>0 and der1[j]>0):
                rc1.append(j)
                r1=len(rc1)
            else:
                r1=0

    r3=0
    if eer3.size!=0:
        for j1 in range(eer3.size):
            j=er3[j1]
            if (der2[i]<0 and der3[j]<0):
                rc3.append(j)
                r3=len(rc3)
            elif (der2[i]>0 and der3[j]>0):
                rc3.append(j)
                r3=len(rc3)
            else:
                r3=0
    r4=0
    if eer4.size!=0:
        for j1 in range(eer4.size):
            j=er4[j1]
            if (der2[i]<0 and der4[j]<0):
                rc4.append(j)
                r4=len(rc4)
            elif (der2[i]>0 and der4[j]>0):
                rc4.append(j)
                r4=len(rc4)
            else:
                r4=0

    if (r1==0 and r3==0) or (r3==0 and r4==0) or (r1==0 and r4==0):
        c2.append(i)

##eliminating events w.r.t. mms3
c3 = []
for i in range(l3):
    ele= smms3[i]
    eer2 = (np.array(np.where(np.logical_and(smms2>=ele-m23, smms2<=ele+m23))))
    eer1 = (np.array(np.where(np.logical_and(smms1>=ele-m13, smms1<=ele+m13))))
    eer4 = (np.array(np.where(np.logical_and(smms4>=ele-m43, smms4<=ele+m43))))
    er2 = list(chain(*eer2))
    er1 = list(chain(*eer1))
    er4 = list(chain(*eer4))
    rc2=[]
    rc1=[]
    rc4=[]
    r2=0
    if eer2.size!=0:
        for j1 in range(eer2.size):
            j=er2[j1]
            if (der3[i]<0 and der2[j]<0):
                rc2.append(j)
                r2=len(rc2)
            elif (der3[i]>0 and der2[j]>0):
                rc2.append(j)
                r2=len(rc2)
            else:
                r2=0

    r1=0
    if eer1.size!=0:
        for j1 in range(eer1.size):
            j=er1[j1]
            if (der3[i]<0 and der1[j]<0):
                rc1.append(j)
                r1=len(rc1)
            elif (der3[i]>0 and der1[j]>0):
                rc1.append(j)
                r1=len(rc1)
            else:
                r1=0
    r4=0
    if eer4.size!=0:
        for j1 in range(eer4.size):
            j=er4[j1]
            if (der3[i]<0 and der4[j]<0):
                rc4.append(j)
                r4=len(rc4)
            elif (der3[i]>0 and der4[j]>0):
                rc4.append(j)
                r4=len(rc4)
            else:
                r4=0

    if (r2==0 and r4==0) or (r4==0 and r1==0) or (r1==0 and r2==0):
        c3.append(i)

##eliminating events w.r.t. mms4
c4 = []
for i in range(l4):
    ele= smms4[i]
    eer2 = (np.array(np.where(np.logical_and(smms2>=ele-m24, smms2<=ele+m24))))
    eer3 = (np.array(np.where(np.logical_and(smms3>=ele-m43, smms3<=ele+m43))))
    eer1 = (np.array(np.where(np.logical_and(smms1>=ele-m14, smms1<=ele+m14))))
    er2 = list(chain(*eer2))
    er3 = list(chain(*eer3))
    er1 = list(chain(*eer1))
    rc2=[]
    rc3=[]
    rc1=[]
    r2=0
    if eer2.size!=0:
        for j1 in range(eer2.size):
            j=er2[j1]
            if (der4[i]<0 and der2[j]<0):
                rc2.append(j)
                r2=len(rc2)
            elif (der4[i]>0 and der2[j]>0):
                rc2.append(j)
                r2=len(rc2)
            else:
                r2=0

    r3=0
    if eer3.size!=0:
        for j1 in range(eer3.size):
            j=er3[j1]
            if (der4[i]<0 and der3[j]<0):
                rc3.append(j)
                r3=len(rc3)
            elif (der4[i]>0 and der3[j]>0):
                rc3.append(j)
                r3=len(rc3)
            else:
                r3=0
    r1=0
    if eer1.size!=0:
        for j1 in range(eer1.size):
            j=er1[j1]
            if (der4[i]<0 and der1[j]<0):
                rc1.append(j)
                r1=len(rc1)
            elif (der4[i]>0 and der1[j]>0):
                rc1.append(j)
                r1=len(rc1)
            else:
                r1=0

    if (r2==0 and r3==0) or (r3==0 and r1==0) or (r2==0 and r1==0):
        c4.append(i)



# indices to eliminate from each mms satellite events 
np.save('2cI1',c1)
np.save('2cI2',c2)
np.save('2cI3',c3)
np.save('2cI4',c4)
