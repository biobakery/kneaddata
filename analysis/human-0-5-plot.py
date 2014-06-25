
# coding: utf-8

# In[1]:

#%matplotlib inline
try:
    import cPickle as pickle
except:
    try:
        import _pickle as pickle
    except:
        import pickle
    
import numpy as np
import matplotlib.pyplot as plt


# In[2]:

with open("human-0-5.npy") as f:
    lstrFnames, lstrFieldNames, aData = pickle.load(f)


# In[3]:

print(lstrFnames)
print(aData)


# In[4]:

index_zipped = zip(lstrFieldNames, range(len(lstrFieldNames)))
print(index_zipped)


# In[5]:

# P(human | it's really human)
aSpecificity = np.divide(aData[:,8] + aData[:,9] + aData[:,10], aData[:,5] + aData[:,6] + aData[:,7])
print(aSpecificity)
aSpecificity[0] = 1
print(aSpecificity)


# In[6]:

# How much of the original data is human reads?
aPercentContam = np.divide(aData[:,4], aData[:,0])
print(aPercentContam)


# In[7]:

# P(nonhuman | nonhuman)
aNonBacterialOrig = (aData[:,1] + aData[:,2] + aData[:,3])-(aData[:,5] + aData[:,6] + aData[:,7])
aNonBacterialRemoved = aData[:,11] + aData[:,12] + aData[:,13]
aSensitivity = np.divide(np.subtract(aNonBacterialOrig, aNonBacterialRemoved),aNonBacterialOrig )
print(aSensitivity)


# In[8]:

# Code for plotting 
plt.figure()
plt.xlim(-0.01,1)
plt.ylim(0.999,1.0001)
plt.plot(aPercentContam, aSpecificity, 'o', color='r', label="P(Human|Human)")
plt.plot(aPercentContam, aSensitivity, '^', color='b', label="P(Nonhuman|Nonhuman)")
plt.title("Specificity and Sensitivity of Pipeline, as a Function of Human Contamination")
plt.xlabel("Contamination (proportion)")
plt.ylabel("Proportion")
plt.legend(loc=4)
plt.show()

