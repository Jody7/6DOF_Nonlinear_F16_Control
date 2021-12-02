local f_counter = 0

local function mf16(alpha, beta, de, da, dr, p, q, r, cbar, b, V, xcg, xcgref)
	f_counter = f_counter + 1
	local phat = p * b / (2 * V)
	local qhat = q * cbar / (2 * V)
	local rhat = r * b / (2 * V)

	local a0 = -1.943367e-2
	local a1 = 2.136104e-1
	local a2 = -2.903457e-1
	local a3 = -3.348641e-3
	local a4 = -2.060504e-1
	local a5 = 6.988016e-1
	local a6 = -9.035381e-1

	local b0 = 4.833383e-1
	local b1 = 8.644627
	local b2 = 1.131098e1
	local b3 = -7.422961e1
	local b4 = 6.075776e1

	local c0 = -1.145916
	local c1 = 6.016057e-2
	local c2 = 1.642479e-1

	local d0 = -1.006733e-1
	local d1 = 8.679799e-1
	local d2 = 4.260586
	local d3 = -6.923267

	local e0 = 8.071648e-1
	local e1 = 1.189633e-1
	local e2 = 4.177702
	local e3 = -9.162236

	local f0 = -1.378278e-1
	local f1 = -4.211369
	local f2 = 4.775187
	local f3 = -1.026225e1
	local f4 = 8.399763
	local f5 = -4.354000e-1

	local g0 = -3.054956e1
	local g1 = -4.132305e1
	local g2 = 3.292788e2
	local g3 = -6.848038e2
	local g4 = 4.080244e2

	local h0 = -1.05853e-1
	local h1 = -5.776677e-1
	local h2 = -1.672435e-2
	local h3 = 1.357256e-1
	local h4 = 2.172952e-1
	local h5 = 3.464156
	local h6 = -2.835451
	local h7 = -1.098104

	local i0 = -4.126806e-1
	local i1 = -1.189974e-1
	local i2 = 1.247721
	local i3 = -7.391132e-1

	local j0 = 6.250437e-2
	local j1 = 6.067723e-1
	local j2 = -1.101964
	local j3 = 9.100087
	local j4 = -1.192672e1

	local k0 = -1.463144e-1
	local k1 = -4.07391e-2
	local k2 = 3.253159e-2
	local k3 = 4.851209e-1
	local k4 = 2.978850e-1
	local k5 = -3.746393e-1
	local k6 = -3.213068e-1

	local l0 = 2.635729e-2
	local l1 = -2.192910e-2
	local l2 = -3.152901e-3
	local l3 = -5.817803e-2
	local l4 = 4.516159e-1
	local l5 = -4.928702e-1
	local l6 = -1.579864e-2

	local m0 = -2.029370e-2
	local m1 = 4.660702e-2
	local m2 = -6.012308e-1
	local m3 = -8.062977e-2
	local m4 = 8.320429e-2
	local m5 = 5.018538e-1
	local m6 = 6.378864e-1
	local m7 = 4.226356e-1

	local n0 = -5.19153
	local n1 = -3.554716
	local n2 = -3.598636e1
	local n3 = 2.247355e2
	local n4 = -4.120991e2
	local n5 = 2.411750e2

	local o0 = 2.993363e-1
	local o1 = 6.594004e-2
	local o2 = -2.003125e-1
	local o3 = -6.233977e-2
	local o4 = -2.107885
	local o5 = 2.141420
	local o6 = 8.476901e-1

	local p0 = 2.677652e-2
	local p1 = -3.298246e-1
	local p2 = 1.926178e-1
	local p3 = 4.013325
	local p4 = -4.404302

	local q0 = -3.698756e-1
	local q1 = -1.167551e-1
	local q2 = -7.641297e-1

	local r0 = -3.348717e-2
	local r1 = 4.276655e-2
	local r2 = 6.573646e-3
	local r3 = 3.535831e-1
	local r4 = -1.373308
	local r5 = 1.237582
	local r6 = 2.302543e-1
	local r7 = -2.512876e-1
	local r8 = 1.588105e-1
	local r9 = -5.199526e-1

	local s0 = -8.115894e-2
	local s1 = -1.156580e-2
	local s2 = 2.514167e-2
	local s3 = 2.038748e-1
	local s4 = -3.337476e-1
	local s5 = 1.004297e-1

	local Cx0 = a0 + a1 * alpha + a2 * math.pow(de, 2) + a3 * de + a4 * alpha * de + a5 * math.pow(alpha, 2) + a6 * math.pow(alpha, 3)
	local Cxq = b0 + b1 * alpha + b2 * math.pow(alpha, 2) + b3 * math.pow(alpha, 3) + b4 * math.pow(alpha, 4)
	local Cy0 = c0 * beta + c1 * da + c2 * dr
	local Cyp = d0 + d1 * alpha + d2 * math.pow(alpha, 2) + d3 * math.pow(alpha, 3)
	local Cyr = e0 + e1 * alpha + e2 * math.pow(alpha, 2) + e3 * math.pow(alpha, 3)
	local Cz0 = (f0 + f1 * alpha + f2 * math.pow(alpha, 2) + f3 * math.pow(alpha, 3) + f4 * math.pow(alpha, 4)) * (1 - math.pow(beta,2)) + f5 * de
	local Czq = g0 + g1 * alpha + g2 * math.pow(alpha, 2) + g3 * math.pow(alpha, 3) + g4 * math.pow(alpha, 4)
	local Cl0 = h0 * beta + h1 * alpha * beta + h2 * math.pow(alpha, 2) * beta + h3 * math.pow(beta,2) + h4 * alpha * math.pow(beta,2) + h5 * math.pow(alpha, 3) * beta + h6 * math.pow(alpha, 4) * beta + h7 * math.pow(alpha, 2) * math.pow(beta,2)
	local Clp = i0 + i1 * alpha + i2 * math.pow(alpha, 2) + i3 * math.pow(alpha, 3)
	local Clr = j0 + j1 * alpha + j2 * math.pow(alpha, 2) + j3 * math.pow(alpha, 3) + j4 * math.pow(alpha, 4)
	local Clda = k0 + k1 * alpha + k2 * beta + k3 * math.pow(alpha, 2) + k4 * alpha * beta + k5 * math.pow(alpha, 2) * beta + k6 * math.pow(alpha, 3)
	local Cldr = l0 + l1 * alpha + l2 * beta + l3 * alpha * beta + l4 * math.pow(alpha, 2) * beta + l5 * math.pow(alpha, 3) * beta + l6 * math.pow(beta,2)
	local Cm0 = m0 + m1 * alpha + m2 * de + m3 * alpha * de + m4 * math.pow(de, 2) + m5 * math.pow(alpha, 2) * de + m6 * math.pow(de, 3) + m7 * alpha * math.pow(de, 2)
	
	local Cmq = n0 + n1 * alpha + n2 * math.pow(alpha, 2) + n3 * math.pow(alpha, 3) + n4 * math.pow(alpha, 4) + n5 * math.pow(alpha, 5)
	local Cn0 = o0 * beta + o1 * alpha * beta + o2 * math.pow(beta, 2) + o3 * alpha * math.pow(beta, 2) + o4 * math.pow(alpha, 2) * beta + o5 * math.pow(alpha, 2) * math.pow(beta, 2) + o6 * math.pow(alpha, 3) * beta
	local Cnp = p0 + p1 * alpha + p2 * math.pow(alpha, 2) + p3 * math.pow(alpha, 3) + p4 * math.pow(alpha, 4)
	local Cnr = q0 + q1 * alpha + q2 * math.pow(alpha, 2)
	local Cnda = r0 + r1 * alpha + r2 * beta + r3 * alpha * beta + r4 * math.pow(alpha, 2) * beta + r5 * math.pow(alpha, 3) * beta + r6 * math.pow(alpha, 2) + r7 * math.pow(alpha, 3) + r8 * math.pow(beta, 3) + r9 * alpha * math.pow(beta, 3)
	local Cndr = s0 + s1 * alpha + s2 * beta + s3 * alpha * beta + s4 * math.pow(alpha, 2) * beta + s5 * math.pow(alpha, 2)
	
	local Cx = Cx0 + Cxq * qhat
	local Cy = Cy0 + Cyp * phat + Cyr * rhat
	local Cz = Cz0 + Czq * qhat
	local Cl = Cl0 + Clp * phat + Clr * rhat + Clda * da + Cldr * dr
	local Cm = Cm0 + Cmq * qhat + Cz * (xcgref - xcg)
	local Cn = Cn0 + Cnp * phat + Cnr * rhat + Cnda * da + Cndr * dr - Cy * (xcgref - xcg) * (cbar / b)
	
	Cx = math.min(math.max(Cx, -0.2), 0.2)
	Cy = math.min(math.max(Cy, -1), 1)
	Cz = math.min(math.max(Cz, -10), 10)
	Cl = math.min(math.max(Cl, -0.2), 0.2)
	Cm = math.min(math.max(Cm, -0.2), 0.2)
	Cn = math.min(math.max(Cn, -0.2), 0.2)
	
	--if f_counter % 30 == 0 then
	--	print(math.round(Cx*1000)/1000,
	--		math.round(Cy*1000)/1000,
	--		math.round(Cz*1000)/1000,
	--		math.round(Cl*1000)/1000,
	--		math.round(Cm*1000)/1000,
	--		math.round(Cn*1000)/1000
	--	)
	--end
	
	return Cx, Cy, Cz, Cl, Cm, Cn
end

return mf16
