"""
========
Barchart
========
A bar plot with errorbars and height labels on individual bars
"""
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import json

SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 20

##plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

algo = 'NSGA3'
function = 'ZDT6'
with open(''.join(['results/HV_',function,'_',algo,'.pk1']),'r') as filename:
    hv_NSGA3 = pickle.load(filename)

Nsamples = 10
NEval = 20000

##NGer = len(hv_NSGA3[0])
##gen_ind = (np.linspace(0,NGer-1,Nsamples+1)[1:]).astype(int)
##eval_ind = (gen_ind + 2)*100

eval_ind = (np.linspace(0,NEval,Nsamples+1)[1:]).astype(int)
gen_ind = eval_ind/100 - 2

meanHV_NSGA3 = hv_NSGA3[:,gen_ind].mean(axis=0)
stdHV_NSGA3 = hv_NSGA3[:,gen_ind].std(axis=0)

ind = np.arange(Nsamples)  # the x locations for the groups
width = 0.2       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind + 3*width, meanHV_NSGA3, width, color='r', yerr=stdHV_NSGA3)

algo = 'NSGA2'
with open(''.join(['results/',algo,'_hvValues_',function,'.json']),'r') as filename:
    hv_NSGA2 = np.asarray(json.load(filename))

meanHV_NSGA2 = hv_NSGA2[:,gen_ind].mean(axis=0)
stdHV_NSGA2 = hv_NSGA2[:,gen_ind].std(axis=0)

rects2 = ax.bar(ind + 0*width, meanHV_NSGA2, width, color='y', yerr=stdHV_NSGA2)

algo = 'SPEA2'
with open(''.join(['results/',algo,'_hvValues_',function,'.json']),'r') as filename:
    hv_SPEA2 = np.asarray(json.load(filename))

meanHV_SPEA2 = hv_SPEA2[:,gen_ind].mean(axis=0)
stdHV_SPEA2 = hv_SPEA2[:,gen_ind].std(axis=0)

rects3 = ax.bar(ind + 1*width, meanHV_SPEA2, width, color='b', yerr=stdHV_SPEA2)

algo = 'SMSEMOA'
with open(''.join(['results/',algo,'_hvValues_',function,'.json']),'r') as filename:
    hv_SMSEMOA = np.asarray(json.load(filename))

meanHV_SMSEMOA = hv_SMSEMOA[:,eval_ind-101].mean(axis=0)
stdHV_SMSEMOA = hv_SMSEMOA[:,eval_ind-101].std(axis=0)

rects4 = ax.bar(ind + 2*width, meanHV_SMSEMOA, width, color='g', yerr=stdHV_SMSEMOA)

# add some text for labels, title and axes ticks
ax.set_ylabel('Average Hypervolume')
ax.set_xlabel('Evaluations')
ax.set_title(''.join(['Average hypervolume evolution for ',function]))
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(eval_ind)
ax.set_ylim([0,max(meanHV_NSGA3)*1.2])

ax.legend((rects2[0],rects3[0],rects4[0],rects1[0]), ('NSGA-II','SPEA2','SMS-EMOA','MOMCEDA'),
           loc='center left', bbox_to_anchor=(1, 0.5),fancybox=True, shadow=True, ncol=1, fontsize='large')

##def autolabel(rects):
##    """
##    Attach a text label above each bar displaying its height
##    """
##    for rect in rects:
##        height = rect.get_height()
##        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
##                '%.2f' % height,
##                ha='center', va='bottom')
##
##autolabel(rects1)
##autolabel(rects2)
##autolabel(rects3)

##plt.show()

plt.savefig(''.join(['HVbar_',function,'.png']), bbox_inches='tight')
