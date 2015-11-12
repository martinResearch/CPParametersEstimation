# Fitting a line paramterized as ax+b=y to a set of points with outliers
# examxple of the method presented in 
# Improving a Constraint Programming Approach for Parameter Estimation
# Bertrand Neveu, Martin de la Gorce and Gilles Trombettoni,
# 27th IEEE International Conference on Tools with Artificial Intelligence 
# Nov 9-11 2015, Vietri sul Mare, Italy 
# Implemented By Martin de La Gorce , novembre 2015
# Ecole des Ponts et Chaussées 

from Qintersection import *

class lineFittingConstraint():
    """this class allows to represent a constraint abs(ax+b-y)<=tolerance 
    tat can be use to find the lines that have fitt a set of observed points x_i,y_i
    with outliers"""
    
    def __init__(self,x,y,tolerance,color):
        self.x=x
        self.y=y
        self.tolerance=tolerance
        self.color=color

    def forwardbackward(self,intervals):    

        if not (isempty(intervals[0]) or isempty(intervals[1])):
            interval_b2=np.sort(self.y-intervals[0]*self.x)
            interval_b2[0]-=self.tolerance
            interval_b2[1]+=self.tolerance   
            intervals[1]=[max(intervals[1][0],interval_b2[0]),
                        min(intervals[1][1],interval_b2[1])]
            
            interval_a2=np.sort((self.y-intervals[1])/self.x)
            interval_a2[0]-=abs(self.tolerance/self.x)
            interval_a2[1]+=abs(self.tolerance/self.x)   
            intervals[0]=np.array([max( intervals[0][0],interval_a2[0]),
                        min( intervals[0][1],interval_a2[1])] ) 
            return intervals

        
    def eval(self,point):
        return abs((self.y-point[0]*self.x)-point[1])<self.tolerance
        
    def slope(self,interval_a):
        return self.x
        
        
    def draw(self,ax,intervals,displayIntervalsBox=False):

        interval_a=intervals[0]
        interval_b=intervals[1]
        if not isempty(interval_a):
            b_1=self.y-interval_a[0]*self.x
            b_2=self.y-interval_a[1]*self.x
            constraintPlot=Polygon([[interval_a[0], b_1-self.tolerance], 
                                    [interval_a[0], b_1+self.tolerance], 
                                    [interval_a[1], b_2+self.tolerance], 
                                    [interval_a[1], b_2-self.tolerance]],
                                   closed=True,color=self.color, alpha=0.3, fill=True)
    
            ax.add_artist(constraintPlot)

            if not isempty(interval_b):
                constraintBox = Rectangle([interval_a[0],interval_b[0]], interval_a[1]-interval_a[0],interval_b[1]-interval_b[0],facecolor="none", edgecolor="none")
                ax.add_artist(constraintBox)
                constraintPlot.set_clip_path(constraintBox)
            if displayIntervalsBox:
                drawBox(ax,intervals,self.color,self.color,alpha=0.1,fill=True)            

        
    def projectConstraintDirection(self,intervals,dir,display=False,ax=None,delta=0):
        interval_a=intervals[0]
        b_1=self.y-interval_a[0]*self.x
        b_2=self.y-interval_a[1]*self.x
        v=np.zeros(4)
        corners=np.zeros((2,4))
        corners[:,0]=[interval_a[0],b_1+self.tolerance]
        corners[:,1]=[interval_a[0],b_1-self.tolerance]
        corners[:,2]=[interval_a[1],b_2+self.tolerance]
        corners[:,3]=[interval_a[1],b_2-self.tolerance]    
        v=dir.dot(corners)
        imin=np.argmin(v)
        imax=np.argmax(v)    
        if display:    
            l1= lines.Line2D([corners[0,imin],1+delta],[corners[1,imin],v[imin]-dir[0]*(1+delta)],color=[0.6,0.6,0.6])
            l2= lines.Line2D([corners[0,imax],1+delta],[corners[1,imax],v[imax]-dir[0]*(1+delta)],color=[0.6,0.6,0.6])
            l3= lines.Line2D([1+delta,1+delta],[v[imin]-dir[0]*(1+delta),v[imax]-dir[0]*(1+delta)],color=[0,0,0])           
            l1.set_clip_on(False)
            l2.set_clip_on(False)
            l3.set_clip_on(False)
            ax.add_line(l1)
            ax.add_line(l2)
            ax.add_line(l3)            
        return np.array([np.min(v),np.max(v)]) 
    
def exampleLineFitting():
       
    nbInliers=3
    Q=nbInliers
    nbOutliers=7
    tolerance=0.01
    epsilon=0.05 # size of the boxes that we consider small enough 
    report=Report('./generated_images/',format='svg')# set to None if you do not want to save images
    #report=None
    np.random.seed(14)
    x_inliers=np.random.rand(nbInliers)
    a=0.3
    b=0.4
    y_inliers=a*x_inliers+b+(np.random.rand(nbInliers)-0.5)*tolerance
   
    x_outliers=np.random.rand(nbOutliers)
    y_outliers=np.random.rand(nbOutliers)
    plt.ion()
   
    plt.xlabel('x')
    plt.xlabel('y')
    plt.axis('tight')
    plt.axis([0,1,0,1])    
   
    x=np.hstack((x_inliers,x_outliers))
    y=np.hstack((y_inliers,y_outliers))
    label=np.hstack((np.ones(nbInliers),np.zeros(nbOutliers)))
    colors='b'*nbInliers+'r'*nbOutliers
    nbPoints=len(x)
   
    fig1=plt.figure(1)
    ax=fig1.gca()
    for i in range(nbPoints):
        circle=plt.Circle((x[i],y[i]),tolerance,color=colors[i])
        ax.add_artist(circle)
    plt.show()
    ax.axis([0,1,0,1])
    axisEqualNoAxesLimitChange(ax)
    if not report is None:
        report.saveFig(ax)
   
    box=np.array([[-1,1],[-1,1]])   

    fig2=plt.figure(2)
    fig2.patch.set_facecolor('#FFFFFF')   
    ax=plt.axes() 
    ax.axis(box.flatten())    
   
    allConstraints=[lineFittingConstraint(x[i],y[i],tolerance,colors[i]) for i in range(nbPoints)]
    qinter=Qintersection(box,allConstraints,Q,epsilon,nbContractions=1,report=report)
    qinter.drawConstraints(ax)
    axisEqualNoAxesLimitChange(ax)
    solutions=qinter.solve(displayBoxes=True,ax=ax) 
    print 'without the projection allong the new direction'
    
    fig3=plt.figure(3)
    fig3.patch.set_facecolor('#FFFFFF')   
    ax=plt.axes() 
    ax.axis(box.flatten())    
    qinter.useProjectionNewDim=False
    solutions=qinter.solve(displayBoxes=True,ax=ax)   
    print 'done'
          

if __name__ == "__main__":
    exampleLineFitting()
       