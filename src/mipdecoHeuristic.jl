export mipdecoHeuristic, heuristicKnapsack

function heuristicKnapsack(mc,gc,delta)
    mt = copy(mc)
    pred = 0
    cIndex = reverse(sortperm(abs.(gc)))
    wCount = 0
	for j=1:length(cIndex)
		ind = cIndex[Int(j)]
			  if  gc[ind] > 0 && mc[ind] == 1
				mt[ind] = 0
				pred -= gc[ind]
                wCount += 1
			  elseif  gc[ind] < 0 && mc[ind] == 0
				mt[ind] = 1
				pred += gc[ind]
                wCount += 1
			  end
            if (wCount >= delta)
                break
            end
	  end
    return mt,pred
end



function  mipdecoHeuristic(mc::BitArray{1}, pInv::InverseParam, pMis;
	                out::Int=2, getNeighborhood::Function=x->trues(size(x)),
					delta::Int=32, sigma::Float64=0.25,outFile::String="")

	mc          = copy(mc)
    iter        = 0
	TermFlag	= 0
	alpha       = pInv.alpha
	maxIter     = pInv.maxIter      #  Max. no. iterations.

	## Time statistics
	timeMIP = 0.0

	## evaluate function and derivatives
    Jc,gc,Dc,Fc,Rc,timeF,timeR = jInvObj(mc,pMis,pInv)
    rId = getNeighborhood(mc)
    gc .*= rId
    J0 = Jc;
	His = zeros(maxIter,3)

	### Initialization for trust region algorithm
	outStr = @sprintf("\n %4s\t%08s\t%08s\t%08s\t%08s\t%08s\t%08s\t%06s\n",
					  	"iter", "Fc", "Rc","alpha","Jc/J0","pred","ared","delta")
    if out>=2; print(outStr); end
	if !isempty(outFile)
		f = open(outFile, "w")
		write(f, outStr)
	end

      outStr = @sprintf("%3d\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3d\n",
		         iter, Fc, Rc,alpha[1],Jc/J0,0,0,delta)
      if out>=2; print(outStr); end
      if !isempty(outFile)
		  write(f, outStr)
	  end

	while (delta > 0)
	  iter  += 1
	  toMIP = @elapsed begin
          mt,pred = heuristicKnapsack(mc,gc,delta)
      end
      pred *=-1
      timeMIP += toMIP

      Jt,gt,Dt,Ft,Rt,toF,toR = jInvObj(mt,pMis,pInv)
      timeF += toF
      timeR += toR
	  ared  = Jc - Jt  	# actual reduction
	  deltaO = delta

      outStr = @sprintf("%3d\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3d\n",
		         iter, Fc, Rc,alpha[1],Jc/J0,pred,ared,delta)
      if out>=2; print(outStr); end
      if !isempty(outFile)
		  write(f, outStr)
	  end

        if ared > 0
	    mc .= mt
	    (Jc,Dc,Fc,Rc) = (Jt,Dt,Ft,Rt)

	    if ared >= sigma*pred;delta *= 2;end
	      gc = gt
          rId = getNeighborhood(mc)
          gc .*= rId
        else
         delta = floor(.5*delta)
	  end
	  His[iter,:] = [Jc Jt deltaO]
	  if iter >= maxIter
          TermFlag = 1
          break
      end
	end 	#while of trust region
	println("Trust region on entire space: MisFit, Regularizer, Total $(Fc) , $(Rc)	, $(Jc)")
	println("Time to compute Data term, Regularizer, and MIP: $(timeF), $(timeR), $(timeMIP) \n")
	if out>=1
		if TermFlag==0
			println("Trust region iterated $iter times with final objective $(Jc)." )
		else
			println("Trust region reached iteration limit with final objective $(Jc).")
		end
	end
    if !isempty(outFile)
	  close(f)
  end
	return mc,Dc,TermFlag,His
end
