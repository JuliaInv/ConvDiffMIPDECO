export mipReg

function mipReg(m::Vector,mref,M::AbstractMesh;Iact=I, q::Int64=1)
	# Rc = e'* m.^3.*(1-m).^3
	Rc = sum(m.^q .* (1.0 .- m).^q)
	dR = q * m.^(q-1) .* (1.0 .- m).^q - q .* m.^q .*(m .- 1.0).^(q-1)
	# d2R = SparseMatrixCSC(Diagonal(max.(0.0 ,- 6.0*m.*(m .- 1.).^3 - 3.0 * m.^3 .* (2.0*m .- 2.) - 18.0*m.^2.0 .* (m .- 1.0).^2)))
	d2R = SparseMatrixCSC(Diagonal( (q*(q-1))*m.^(max(0,q-2)) .*(1.0 .- m).^q -
	                                (q*(q-1)) * m.^q .* (m .- 1.0).^(max(0,q-2))  -
								     2*(q*q)*m.^(q-1) .* (m .- 1.0).^(q-1)))
	d2R = max.(d2R,0.0)
	return Rc,dR,d2R
end
