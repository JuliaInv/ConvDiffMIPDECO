export naiveRounding, objGapRedRounding, massPreservingRounding

"""
naiveRounding(s,thresh=0.5)

simple rounding scheme based on threshold.

"""
naiveRounding(s,thresh=0.5) = (s.>=thresh)

"""
function massPreservingRounding(s)

mass-preserving rounding

"""
function massPreservingRounding(s)
    r   = falses(length(s))
    pos = reverse(sortperm(s))
    for k=1:round(Int64,sum(s))
        r[pos[k]] = 1
    end
    return r
end

"""
function objGapRedRounding(s,pMis,pInv;delta=0.1)

objective gap reduction rounding
"""
function objGapRedRounding(s,pMis,pInv;delta=0.1,out=0)
    r     = zeros(Int64,length(s))
    Jbest = 10000
    sbest = 0
    smin = floor(minimum(s)*10)*.1 + delta
    smax = floor(maximum(s)*10)*.1
    for k=smin:delta:smax
        r = 1.0*naiveRounding(s,k)
        sig,dsig = pInv.modelfun(r)
         Dc,F,dF,d2F,pMis,tMis = computeMisfit(sig,pMis,true)
         R,dR,d2R = computeRegularizer(pInv.regularizer,1.0*r,pInv.mref,pInv.MInv,pInv.alpha)
		 if out>=2
			 println("Obj-gap reduction rounding: cut-off, MisFit, Regularizer, Total $(k), $(F) , $(R) , $(F+R)")
		 end
		 
         if Jbest > F+R
             sbest = k
             Jbest = F+R
         end
     end
     if sbest > 0
         r = naiveRounding(s,sbest)
     end
    return r
end