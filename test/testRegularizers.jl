using ConvDiffMIPDECO
using jInv.InverseSolve
using Test
using jInv.Mesh
using jInv.Utils
using LinearAlgebra

# build regular mesh and Iact
domain = [0;5;0;5;]
n     = [12;13]
M     = getRegularMesh(domain,n)
mc    = randn(M.nc)

regFuns = [(m,mref,M)->mipReg(m,mref,M,q=1), (m,mref,M)->mipReg(m,mref,M,q=3)]
for k=1:length(regFuns)
	# checkDerivative of $(regFuns[k])
	function testFun(x,v=[])
		Sc,dS,d2S = regFuns[k](x,0.0.*x,M)
		if isempty(v)
			return Sc
		else
			return Sc,dot(dS,v)
		end
	end
	chkDer, = checkDerivative(testFun,mc,out=true)
	@test chkDer
end
