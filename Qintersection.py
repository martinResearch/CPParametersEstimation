# Python Implementation of the method presented in 
# Improving a Constraint Programming Approach for Parameter Estimation
# Bertrand Neveu, Martin de la Gorce and Gilles Trombettoni,
# 27th IEEE International Conference on Tools with Artificial Intelligence 
# Nov 9-11 2015, Vietri sul Mare, Italy 
# Implemented By Martin de La Gorce , novembre 2015
# Ecole des Ponts et Chaussees 


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import  Polygon,Rectangle
from matplotlib import pyplot, lines,transforms
from copy import copy,deepcopy
import time
import os


class Report():
    """small class to handle export to images with automatic increment of the image id"""
    def __init__(self,folder,format='svg'):
        self.idimage=0
        self.folder=folder
        if not os.path.isdir(folder):
            os.mkdir(folder)
        self.format=format
    def saveFig(self,ax):
       
        ax.figure.savefig(os.path.join(self.folder,'image%d.%s'%(self.idimage,self.format)))        
        self.idimage+=1        



class QintersectionNode():
    """This class implements a node inthe tree used to explore the parmter space"""
    def __init__(self,allConstraints,box,possibleCS,colors=[],report=None):
        """allConstraints: the list of all observed data constraints
           box: the box in the peramters space 
           possibleCS: the list of observations indices that are kept as possible element of the consensus set"""
        self.box=np.array(box)
        self.allConstraints=allConstraints
        self.colors=colors
        self.possibleCS=copy(possibleCS)
        self.validCS=[]
        self.ndim=len(self.box)
        nbPoints=len(possibleCS)
        self.intervals=np.zeros((nbPoints,self.ndim,2)) 
        self.intervals[:,:,:]=self.box[None,:,:]
        self.report=report
        self.bestPoint=None
            
    def drawConstraints(self,ax):
        for i,idc in enumerate(self.possibleCS):
            self.allConstraints[idc].draw(ax,self.intervals[i],displayIntervalsBox=True)  
        
    def  drawConstraintAndIntervals(self,ax):        
        ax.clear()      
        self.drawConstraints(ax)
        drawProjectedIntervals(ax,self.intervals) 
        self.draw(ax)            
        plt.draw()
        plt.show()
        axisEqualNoAxesLimitChange(ax)
        plt.show()
        if not self.report is None:
            self.report.saveFig(ax)
        return ax   
        
    def contractAndIntersect(self,Q,useProjectionNewDim=True,display=False,ax=None): 
        for i,idconstraint in enumerate(self.possibleCS):   
            self.intervals[i]=self.allConstraints[idconstraint].forwardbackward(self.intervals[i]) 
        if display :
            self.drawConstraintAndIntervals(ax)        
        qbox=np.zeros(self.box.shape)
        for idim in range(self.ndim):
            qbox[idim],self.intervals[:,idim] = Qintersection1D(self.intervals[:,idim],Q)
            if display :
                self.drawConstraintAndIntervals(ax)            
        self.box=qbox 
        
        keep=np.all(self.intervals[:,:,1]>=self.intervals[:,:,0],axis=1) 
        
        if useProjectionNewDim:
            # use projection along new direction
            slopes=[]   
            for i,idc in enumerate(self.possibleCS):
                slopes.append(self.allConstraints[idc].slope(self.intervals[i]))
            dir=np.array([np.mean(slopes),1])
            nbPoints=len(self.possibleCS)        
            intervals_projected=np.zeros((nbPoints,2))
            epsilon=0.01*(self.box[0,1]-self.box[0,0])
            if display:
                ax.clear()
                self.drawConstraints(ax)
                axisEqualNoAxesLimitChange(ax)
            for i,idc in enumerate(self.possibleCS):
                    delta= (i+2)*epsilon   
                    intervals_projected[i]=self.allConstraints[idc].projectConstraintDirection(self.intervals[i],dir,display,ax=ax,delta=delta)              
            qinterval,intervals_projected2 = Qintersection1D(intervals_projected,Q) 
            keep=keep & (intervals_projected2[:,1]>=intervals_projected2[:,0])
            
            if display and self.ndim==2:
                
                b1=(qinterval[0]-qbox[0,0]*dir[0])/dir[1]
                b2=(qinterval[1]-qbox[0,0]*dir[0])/dir[1]
                b3=(qinterval[1]-qbox[0,1]*dir[0])/dir[1]
                b4=(qinterval[0]-qbox[0,1]*dir[0])/dir[1]
                ax.add_patch(Polygon([[qbox[0,0], b1],  [qbox[0,0], b2], [qbox[0,1], b3],  [qbox[0,1], b4]],
                                     closed=True, color='g', alpha=0.5,fill=True))
                plt.draw()
                plt.show()            
                axisEqualNoAxesLimitChange(ax)
                plt.show() 
                if not self.report is None:
                    self.report.saveFig(ax)
        
        #update possibleCS
        new_possibleCS=[]
        for i,p  in enumerate(self.possibleCS):
            if keep[i] :
                new_possibleCS.append(p)
        self.possibleCS=new_possibleCS
    
    def isempty(self):
        return len(self.possibleCS)==0
        
    def isSolution(self):
        return len(self.possibleCS)==len(self.validCS)
        
    def draw(self,ax,color='k'): 
        drawBox(ax,self.box,edgecolor=color,facecolor=None,alpha=1,fill=False)
        plt.show() 
        
    def bisect(self):
        """bisect along the dimension with largest width"""
        width=np.zeros(len(self.box))
        for i in range(len(self.box)):
            width[i]=self.box[i][1]-self.box[i][0]
        dimension_split=np.argmax(width)
        mid=0.5*(self.box[dimension_split][0]+self.box[dimension_split][1])
        
        box1=deepcopy(self.box)
        box1[dimension_split]=[box1[dimension_split][0],mid]        
        node1=QintersectionNode(box=box1,allConstraints=self.allConstraints,possibleCS=self.possibleCS)
        if np.all(self.bestPoint>=box1[:,0])  and np.all(self.bestPoint<=box1[:,1]):
            node1.bestPoint=self.bestPoint
            node1.validCS=self.validCS
        box2=deepcopy(self.box)
        box2[dimension_split]=[mid,box2[dimension_split][1]]
        node2=QintersectionNode(box=box2,allConstraints=self.allConstraints,possibleCS=self.possibleCS)
        if np.all(self.bestPoint>=box2[:,0])  and np.all(self.bestPoint<=box2[:,1]):
            node2.bestPoint=self.bestPoint
            node2.validCS=self.validCS
        return node1,node2
    
    def center(self):
        return np.mean(self.box,axis=1)
    
    def validate(self):
        new_validCS=[]
        center=self.center()
        for idc in self.possibleCS:
            if self.allConstraints[idc].eval(center):
                new_validCS.append(idc)
        if len( new_validCS)>=len(self.validCS):
            self.bestPoint=center
            self.validCS=new_validCS
                
    def width(self):
        width=np.zeros(len(self.box))
        for i in range(len(self.box)):
            width[i]=self.box[i][1]-self.box[i][0] 
        return np.max(width)
              
class Qintersection():
    
    """This class implements the method presented in 
    Improving a Constraint Programming Approach for Parameter Estimation
    Bertrand Neveu, Martin de la Gorce and Gilles Trombettoni,
    27th IEEE International Conference on Tools with Artificial Intelligence 
    Nov 9-11 2015, Vietri sul Mare, Italy """
    
    def __init__(self,box,allConstraints,Q,epsilon,nbContractions=1,useProjectionNewDim=True,report=None): 
        self.allConstraints=allConstraints
          
        self.Q=Q
        self.epsilon=epsilon
        self.box=box
        self.nbContractions=nbContractions
        self.useProjectionNewDim=useProjectionNewDim
        self.report=report
        
    def iterate(self): 
        global idimage
        node=self.nodes.pop() 
        # contract the node        
        display=(self.nbiter==0  and self.displayBoxes)  
        box=node.box
        for k in range(self.nbContractions):
            node.contractAndIntersect(self.Q,display=display,ax=self.ax,useProjectionNewDim=self.useProjectionNewDim)
        if self.nbiter==0 and self.displayBoxes:
            self.ax.clear()
            self.drawConstraints(self.ax)            
        if  self.displayBoxes:
            node.draw(self.ax)
            plt.show()
            if not self.report is None:
                self.report.saveFig(self.ax)
             
        if not node.isempty():
            node.validate()
            if node.width()<self.epsilon or node.isSolution():
                self.solutions.append(node)
                if self.displayBoxes:
                    p=drawBox(self.ax,node.box,edgecolor='k',facecolor='g',alpha=0.8,fill=True)
                    p.set_zorder(3)
                    if not self.report is None:
                        self.report.saveFig(self.ax)                   
            else:
                node1,node2=node.bisect()
                self.nodes.append(node1)
                self.nodes.append(node2)
                if  self.displayBoxes:
                    node1.draw(self.ax)
                    node2.draw(self.ax)
                    plt.show()
                    if not self.report is None:
                        self.report.saveFig(self.ax)
        else:
                p=drawBox(self.ax,box,edgecolor='k',facecolor='k',alpha=0.3,fill=True)
                p.set_zorder(3)
                if not self.report is None:
                    self.report.saveFig(self.ax)  
        if self.displayBoxes:
            time.sleep(0.5)  
                
    def solve(self,displayBoxes=False,ax=None):
        self.solutions=[]
        node=QintersectionNode(box=self.box,allConstraints=self.allConstraints,possibleCS=np.arange(len(self.allConstraints)))
        self.nodes=[node]         
        self.displayBoxes=displayBoxes
        self.ax=ax
        self.nbiter=0
        while len(self.nodes)>0:
            self.iterate()
            self.nbiter+=1
        print 'found %d solutions after %d iterations'%(len(self.solutions),self.nbiter)
        return self.solutions
    
    def drawConstraints(self,ax):        
        for constraint in self.allConstraints:
            constraint.draw(ax,self.box)        

def axisEqualNoAxesLimitChange(ax):
    """plt.axis('equal') changes the limits of you axes, 
    this fonction set the aspect ration in pixel to 1 by 
    changing the position of the axis in the figure keeping the limits unchanged
    Martin de La Gorce 
    """
    ax.set_position([0.1,0.1,0.8,0.8])
    t= ax.transData.transform([[0,0],[1,0],[0,1]])
    r=(t[2,1]-t[0,1])/(t[1,0]-t[0,0])
    if r<1:
        width=0.8*r
        pos2 = [(1-width)/2, 0.1,  width, 0.8] 
    else:
        pos2 = [0.1, 0.1,  0.8, 0.8/r]
    ax.set_position(pos2) 
    plt.draw()
    plt.show()

def isempty(interval):
    """we represent empty intervals by intervals whose upper bound is lower than the lower bound
    this make vectorization of the code using numpy easier than handling list of intervals with 
    potentialy empty intervals represented by None"""
    return interval[1]<interval[0]


def drawBox(ax,box,edgecolor,facecolor,alpha,fill):
    """Draws a rectangulare box"""
    if box.shape[1]!=2:
        print 'not yet coded with dimension different from two'
        raise
    if not (isempty(box[0]) or isempty(box[1])):
        patch=Polygon([[box[0][0], box[1][0]],[box[0][0], box[1][1]], 
                              [box[0][1], box[1][1]],[box[0][1], box[1][0]]],
                             closed=True,edgecolor=edgecolor, facecolor=facecolor, alpha=alpha, fill=fill)
        ax.add_patch(patch)
    else :
        patch=None
    
    return patch
       



def drawProjectedIntervals(ax,intervals):
    """this function draw line segment on the left and the top of the image 
    where each segment represents the interval corresponding to 
    and interval that approximates the projection of a constraints 
    on one of the two dimension obtained by forward-backward 
    this 1D segments are then used in the 1 dimensional q intersection function"""
    if intervals.shape[1]!=2:
        print 'not coded when the dimension of the paramter space is different from 2'
        raise
    intervals_a=intervals[:,0]
    axis=ax.axis()
    epsilon=0.01*(axis[3]-axis[2])
    for i in range(len(intervals_a)):
        if not isempty(intervals_a[i]) :
            line = lines.Line2D(intervals_a[i],[axis[3]+(i+2)*epsilon]*2,color='k')
            line.set_clip_on(False)
            ax.add_line(line)      
    intervals_b=intervals[:,1]
    for i in range(len(intervals_b)):
        if not intervals_b[i] is None:
            line = lines.Line2D([axis[1]+(i+2)*epsilon]*2,intervals_b[i],color='k')
            line.set_clip_on(False)#make the line outside the axis still visible
            ax.add_line(line)     
    
 
        
def Qintersection1D(intervals,Q,display=False):
    """This function implement the Q intersection operator one 1-dimensional intervals,
    which can be done exaplaty in nlog(n) with n the number of intervals"""
    onsets=[]
    for interval in intervals:
        if not isempty(interval):      
            onsets.append((interval[0],1))
            onsets.append((interval[1],-1))
    
    onsets=np.array(onsets)
    if len(onsets)>0:
        order=np.argsort(onsets[:,0])
        onsets=onsets[order,:]
        c=np.cumsum(onsets[:,1])
        v=np.nonzero(c>=Q)[0]
        assert(c[-1]==0)
    else :
        v=[]
    if len(v)>0:
        start=onsets[v[0],0]
        end=onsets[v[-1]+1,0]
    else:
        start=np.inf
        end=-np.inf
    
    if display:
        fig=plt.figure()
        ax2=fig.gca()
        c2=np.vstack((np.hstack((0,c)),np.hstack((0,c)))).T.flatten()
        ax2.plot(np.tile(onsets[:,0],(2,1)).T.flatten(),c2[1:-1],linewidth=3)
        h=np.max(c)+1
        ax2.axis([-1,1,0,h])
        drawProjectedIntervals_a(ax2,intervals_a) 
        ax2.plot([-1,1],[Q,Q],':g',linewidth=2)
        ax2.plot([start,start],[0,h],'r:',linewidth=2)
        ax2.plot([end,end],[0,h],'r:',linewidth=2)
        global idimage
        plt.show()
        fig.savefig('image%d.png'%idimage)
        idimage+=1 
    
    new_intervals=np.zeros(intervals.shape)
    for i,interval in enumerate(intervals):
        if  not isempty( interval) :
            new_intervals[i]=np.array([max(interval[0],start),min(interval[1],end)])   
        else:
            
            new_intervals[i]=interval
    return np.array([start,end]),new_intervals

