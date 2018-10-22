import numpy as np

def getSimNames(folderName, Nbins, maxQ, Iebeam, screenDist, elEnergy):
  fileNameSuffix = "_Qmax-{:.6f}".format(maxQ)\
      + "_Ieb-{:.6f}".format(Iebeam)\
      + "_scrnD-{:.6f}".format(screenDist)\
      + "_elE-{:.6f}".format(elEnergy)\
      + "_Bins[" + str(Nbins) + "].dat";

  atmDiffFile = folderName + "/nitrobenzene_atmDiffractionPatternLineOut"\
                  + fileNameSuffix
  molDiffFile = folderName + "/nitrobenzene_molDiffractionPatternLineOut"\
                  + fileNameSuffix

  return atmDiffFile, molDiffFile



def getZeroXQ(Q, molDiff):

  search = molDiff[:-1]*molDiff[1:]
  inds = np.where(search <= 0)
  return Q[inds]


def findzeros(func,ind):
    Y = np.sign(func)
    x0 = []
    #print func[5000]
    #for i in np.arange(len(Y)):
    #    if Y[i]==0:
    #        x0.append(i)
    #Y[x0] = 1
    n0 = []
    for i in np.arange(len(ind)-1):
        if ~np.isnan(Y[i+1]+Y[i]) and abs(Y[i+1]-Y[i])>=2:
            n0.append(i)
            #print Y[i+1], Y[i], i, Y[4998:5001],func[4998:5001]
    #print n0
    for i in np.arange(len(n0)):
        if abs(func[n0[i]])>abs(func[n0[i]+1]):
            n0[i] = n0[i] + 1
    return n0


def fit_scatt_bkg(null,nullst,Q):
    
    # Make sure all null are >0
    null = null[null>0]
    nullst = nullst[null>0]
    
    iters = len(null)-2
    
    # Initialize arrays:
    A = np.zeros((iters,))
    B = np.zeros_like(A)
    C = np.zeros_like(A)
    fitValues = np.zeros((len(nullst),iters))
    BkgFits = np.zeros((len(Q),iters))
    BkgCurve = np.zeros_like(Q)
    
    # Loop over different fit sections
    for i in np.arange(iters):
        # create parameters for function f:
        s1 = nullst[i]
        s2 = nullst[i+1]
        s3 = nullst[i+2]
        alpha1 = np.log(null[i])
        alpha2 = np.log(null[i+1])
        alpha3 = np.log(null[i+2])
        
        # Initialize variables for finding nullstellen in f
        czeroindex=[]
        zz=2
        cstep = 0.00001/zz
        cstart = 0.05/zz
        

        while czeroindex==[] and cstart<500:
            cstep=cstep*zz
            #print cstep
            cstart=cstart*zz
            c = np.arange(-1*cstart,cstart,cstep)
            f = np.power(s3,c)*(alpha2-alpha1)-0.5*(alpha2-alpha1)*(np.power(s1,c)+np.power(s2,c))- \
                (np.power(s2,c)-np.power(s1,c))*(alpha3-alpha1/2-alpha2/2)
            f[np.where(np.round(f,7)==0)] = 0

            cindeces = findzeros(f,c)
            #print cindeces
            if cindeces != []:
                czeros = (np.array(cindeces)-1)*cstep-cstart
                czeroindex = np.where(abs(czeros) > cstep)[0]
                #print czeroindex
            else:
                czeros = []
           
        print "here"
        print czeroindex
        print "here"
        if czeroindex.size == 1:
            C[i] = czeros[czeroindex]
            B[i] = (alpha2-alpha1)/(np.power(s2,C[i])-np.power(s1,C[i]))
            A[i] = 0.5*(alpha1 + alpha2 -B[i]*(np.power(s1,C[i])+np.power(s2,C[i])))
            fitValues[:,i] = np.exp(A[i]+B[i]*np.power(nullst,C[i]))
            BkgFits[:,i] = np.exp(A[i]+B[i]*np.power(Q,C[i]))
            
            #print czeros[czeroindex],C[i],B[i],A[i]
            #print BkgFits[:,i] 
        else:
            BkgFits[:,i] = np.zeros((1,len(Q)))
        
    sRange = np.where(Q<=nullst[1])
    #print sRange
    if ('A' in locals() or 'A' in globals()):
        if len(np.where(BkgFits[:,0]==0)[0])==0:
            BkgCurve[sRange] = BkgFits[sRange,0]
            #print BkgCurve[sRange]
    
            if iters>1:
                for j in np.arange(1,iters):
                    sRange = np.where((Q<=nullst[j+1])&(Q>nullst[j]))
                    BkgCurve[sRange] = (nullst[j+1]-Q[sRange])/(nullst[j+1]-nullst[j])*BkgFits[sRange,j-1] + \
                                       (Q[sRange]-nullst[j])/(nullst[j+1]-nullst[j])*BkgFits[sRange,j]
                    #print BkgCurve[sRange]

                sRange = np.where(Q>nullst[j])
                BkgCurve[sRange] = BkgFits[sRange,j]
            else:
                BkgCurve[:len(Q)] = BkgFits[:,0]
        else:
            BkgCurve[sRange] = 0
    else:
        BkgCurve[sRange] = 0
    
    if len(BkgCurve)!=len(Q):
        BkgCurve = 0
    return BkgCurve


