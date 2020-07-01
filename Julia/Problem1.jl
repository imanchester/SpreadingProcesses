using Mosek
using MAT
using SparseArrays


printstream(msg::AbstractString) = print(msg)

#set parameters
n=1000; 
delta=0.5;
Gamma=200;
dr=3.1;
betaL=1e-8;

#import data needed (Beta matrix, Cost, Likelihood, Vegetation, adjacency matrix) 
file3 = matopen("P1JuliaNew.mat") 
Beta=read(file3,"Beta") #transition matrix
C=read(file3,"C") #cost
Lambda=read(file3,"Lambda") #Likelihood
Pveg=read(file3,"Pveg") #Vegetation
AM = read(file3,"AM") #adapted adjacency matrix (Atemp in Matlab files)
close(file3)

(II,JJ,VV)=findnz(AM)
(LL,MM,WW)=findnz(Beta)
nR= length(VV) #number of nonzero r_ij elements  
nR_slack=length(WW) #number of nonzero Beta elements


#order of variables [slack; y_i; r_ij; exp Cone variables; 
#basic dimensions
numvar=n+nR+1; #+1 for max objective 
numLinCon=1+n;  #+n for log(x(0) <= slack - y <= Inf)
numExp= nR_slack+2*n #(no of i=8 +2)*j 

#linear part of the problem
cval=zeros(numvar)
cval[1]=1
asubi=repeat(1:n, inner=2)    #row index
asubi=append!(asubi,ones(nR)*(n+1))# add sum rij <= Gamma
asubj=ones(2*n)   
asubj[2:2:2*n]=2:n+1 #column index 
asubj=append!(asubj,n+2:n+1+nR)
aval=repeat([1,-1],outer=n)
aval=append!(aval,ones(nR))    #values
ajval=[1, 1, 1, -1]
ajvalN=[1, 1, -1]
aj2val=ones(n)
aj3val=ones(2*n)

q=25 #row length 

dstop=numvar+3*n 
cstop=dstop+3*n 

a2i=Array(numLinCon+1:numLinCon+n)
a2c=Array(numvar+3:3:dstop)

a3i=repeat(numLinCon+n+1:numLinCon+2*n,inner=2) 
a3c=repeat(2:n+1,inner=2)
a3c[1:2:2*n]= Array(dstop+3:3:cstop) #add slack column variables 


Loopstartrow=numLinCon+2*n+1
Loopstartcol=cstop+3 
r_pos=n+2
#setting up linear part of logsumexp

z=Loopstartcol
c = numLinCon + numExp + 1
expStart = numvar + 1

RIJ=Loopstartrow*ones(4)


maketask() do task
	putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))
	global z, r_pos, Loopstartrow, c
	appendcons(task,numLinCon + numExp + n) # Append 'numcon' empty constraints.
    appendvars(task,numvar + 3*numExp) # Append 'numvar' variables. (per cone (ui,1=ti,zi))

    #set up objective 
    #putcslice(task,1, numvar, cval)
    putcj(task,1,1)
    putvarboundsliceconst(task,1, n+2, MSK_BK_FR , -Inf , +Inf ) #no constraints on slack1, y
    putvarboundslice(task,n+2,numvar+1,[MSK_BK_RA for i in 1:nR], [0 for i in 1:nR], [log(VV[i]/betaL) for i in 1:nR]) #add constraints on rij
      

    #set up linear constraints
   	putaijlist(task,asubi, asubj, aval)
    putconboundslice(task,1, numLinCon, [MSK_BK_LO for i in 1:n], [log(Lambda[i]) for i in 1:n], [+Inf for i in 1:n]) 
    putconbound(task,numLinCon, MSK_BK_RA, 0, Gamma)

    # set up linear constraints in exp() 
    putaijlist(task,a2i,a2c,aj2val)
    putconboundsliceconst(task,numLinCon+1,numLinCon+n+1,MSK_BK_FX ,log((1-delta)/(1+dr)) ,log((1-delta)/(1+dr)) ) #j times

    putaijlist(task,a3i,a3c,aj3val) 
    putconboundslice(task,numLinCon+n+1,numLinCon+2*n+1,[MSK_BK_FX for i in 1:n],[log(C[i]/(1+dr)) for i in 1:n],[log(C[i]/(1+dr)) for i in 1:n])
    
    RNTeller=cstop+1

    for j=1:n
    	#determine i values
    	teller=8
    	ii=[j-1-q,j-q,j+1-q,j-1,j+1,j-1+q,j+q,j+1+q] #q is number of rows

    	if j < q+1 
    		ii[1:3].=0#set first 3 values to zero
    	end
    	if j > n-q
    		ii[6:8].=0#set last 3 values to zero
    	end
    	if round(j/q)==j/q  #bottom row
    		ii[[3,5,8]].=0#remove values 3, 5, 8
    	end
    	if round((j+q-1)/q)== (j+q-1)/q #top row
    		ii[[1,4,6]].=0#remove values 1, 4, 6
	   	end
	   	
	   	for i in ii 
	   		if i >0 && Beta[i,j]>0
	   			if Pveg[j]>0
					putaijlist(task,Loopstartrow*ones(4), [z,r_pos, j+1, i+1],ajval) 
					r_pos=r_pos+1 	 
     			else
     				putaijlist(task,Loopstartrow*ones(3), [z, j+1, i+1],ajvalN)
     			end
     			z=z+3
     			putconbound(task,Loopstartrow,MSK_BK_FX,log(Beta[i,j]/(1+dr)),log(Beta[i,j]/(1+dr))) #for each i, for each j   
     			Loopstartrow=Loopstartrow+1
     		else 
     			teller=teller-1
     		end
    	end

    	if teller>0
	    	#set up u as constraint
	    	putarow(task,c, [expStart+3*(j-1); dstop+1+3*(j-1); Array(range(RNTeller, step=3, stop=RNTeller+3*teller-1))], ones(teller+2)) 
	    	putconbound(task,c, MSK_BK_FX, 1.0, 1.0) 
	    	c=c+1
	    	RNTeller=RNTeller+3*teller
    	end

    end 


    #add sum u as constraint
    #c = numLinCon + numExp + 1
    #expStart = numvar + 1
	putvarboundlistconst(task,Array(range(expStart, step=3, stop=expStart + 3*numExp-1)), MSK_BK_FR , -Inf , +Inf ) #no bound on u

    # z_i unbounded
	putvarboundlistconst(task,Array(range(expStart + 2, step=3, stop=expStart + 1 + 3*numExp)), MSK_BK_FR , -Inf , +Inf) 
    # t_i = 1
    #putvarboundlist(task,Array(range(expStart + 1, step=3, stop=expStart + 3*numExp)), [MSK_BK_FX for i in 1:numExp], ones(numExp), ones(numExp))
    putvarboundlistconst(task,Array(range(expStart + 1, step=3, stop=expStart + 3*numExp)), MSK_BK_FX , 1.0, 1.0) 

    #append cones; Every triple is in an exponential cone
    appendconesseq(task,[MSK_CT_PEXP for i in 1:numExp], zeros(numExp), 3*ones(numExp), expStart) 

    #input objective (min/max)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)

    #optimize task
	optimize(task)
	putintparam(task,MSK_IPAR_INFEAS_REPORT_AUTO, MSK_ON)
	solutionsummary(task,MSK_STREAM_LOG)


	prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)
    x = getxx(task,MSK_SOL_ITR);

    if solsta == MSK_SOL_STA_OPTIMAL
        println("Optimal solution")
		#solution summary
		x = getxx(task,MSK_SOL_ITR);
		KK= x[n+2:n+1+nR];
		display(sum(KK))
		#export solution to mat file
		file2 = matopen("RA_N.mat", "w")
		write(file2, "KK", KK)
		close(file2)
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        println("Primal or dual infeasibility.")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        println("Primal or dual infeasibility.")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        println("Other solution status")
    end


	

end






