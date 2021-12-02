local morellif16 = require("morellif16.lua")

local function t_copy(original)
	local copy = {}
	for key, value in pairs(original) do
		copy[key] = value
	end
	return copy
end

local function adc(vt, alt)
	--get mach num and dyn pressure
	local ro = 2.377e-3
	local tfac = 1 - .703e-5 * alt
	local t = 0

	if alt >= 35000 then
		t = 390
	else
		t = 519 * tfac
	end

	local rho = ro * math.pow(tfac, 4.14)
	local a = math.sqrt(1.4 * 1716.3 * t)
	local amach = vt / a
	local qbar = .5 * rho * vt * vt
  
	return amach, qbar
end
	
local function tgear(thtl)
	local tg
	if (thtl <= .77) then
		tg = 64.94 * thtl
	else
		tg = 217.38 * thtl - 117.38
	end
	return tg
end
	
local function rtau(dp)
	local rt
	if (dp <= 25) then
		rt = 1.0
	elseif (dp >= 50) then
		rt = 0.1
	else
		rt = 1.9 - .036 * dp
	end
	return rt
end
	
local function pdot(p3, p1)
		local p2
		if (p1 >= 50) then
			if (p3 >= 50) then
				t = 5
				p2 = p1
			else
				p2 = 60
				t = rtau(p2 - p3)
			end
		else
			if (p3 >= 50) then
				t = 5
				p2 = 40
			else
				p2 = p1
				t = rtau(p2 - p3)
			end
		end
		local pd = t * (p2 - p3)
		return pd
end	
	
local function fix(ele)
	if ele > 0 then
		return math.floor(ele)
	else
		return math.ceil(ele)
	end
end	

local function thrust(power, alt, rmach)
	local a = 
		{{ 1060,   635,    60, -1020, -2700, -3600,},
	{  670,   425,    25,  -170, -1900, -1400,},
	{  880,   690,   345,  -300, -1300,  -595,},
	{ 1140,  1010,   755,   350,  -247,  -342,},
	{ 1500,  1330,  1130,   910,   600,  -200,},
	{ 1860,  1700,  1525,  1360,  1100,   700,}}

	local b = 
		{{12680, 12680, 12610, 12640, 12390, 11680,},
	{ 9150,  9150,  9312,  9839, 10176,  9848,},
	{ 6200,  6313,  6610,  7090,  7750,  8050,},
	{ 3950,  4040,  4290,  4660,  5320,  6100,},
	{ 2450,  2470,  2600,  2840,  3250,  3800,},
	{ 1400,  1400,  1560,  1660,  1930,  2310,}}

	local c = 
		{{20000, 21420, 22700, 24240, 26070, 28886,},
	{15000, 15700, 16860, 18910, 21075, 23319,},
	{10800, 11225, 12250, 13760, 15975, 18300,},
	{ 7000,  7323,  8154,  9285, 11115, 13484,},
	{ 4000,  4435,  5000,  5700,  6860,  8642,},
	{ 2500,  2600,  2835,  3215,  3950,  5057,}}

	if alt < 0 then
		alt = 0.01
	end
	local h = .0001 * alt
	local i = fix(h)
	if (i >= 5) then
		i = 4
	end
	local dh = h- i
	local rm = 5 * rmach
	local m = fix(rm)

	if (m >= 5) then
		m = 4
	elseif (m <= 0) then
		m = 0
	end

	local dm = rm - m
	local cdh = 1 - dh

	i = i + 1
	m = m + 1

	local s = b[i][m] * cdh + b[i + 1][m] * dh
	local t = b[i][m + 1] * cdh + b[i + 1][m + 1] * dh
	local tmil = s + (t - s) * dm

	local tidl = 0
	local thrst = 0
	local tmax = 0

	if (power < 50) then
		s = a[i][m] * cdh + a[i + 1][m] * dh
		t = a[i][m + 1] * cdh + a[i + 1][m + 1] * dh
		tidl = s + (t - s) * dm
		thrst = tidl + (tmil - tidl) * power * .02
	else
		s = c[i][m] * cdh + c[i + 1][m] * dh
		t = c[i][m + 1] * cdh + c[i + 1][m + 1] * dh
		tmax = s + (t - s) * dm
		thrst = tmil + (tmax - tmil) * (power - 50) * .02
	end

	return thrst
end	

local function dampp(alpha)
	local a =
		{{ -0.267,   0.882,  -0.108,  -8.8,    -0.126,  -0.36,   -7.21,   -0.38,    0.061},
	{ -0.11,    0.852,  -0.108, -25.8,    -0.026,  -0.359,  -0.54,   -0.363,   0.052},
	{  0.308,   0.876,  -0.188, -28.9,     0.063,  -0.443,  -5.23,   -0.378,   0.052},
	{  1.34,    0.958,   0.11,  -31.4,     0.113,  -0.42,   -5.26,   -0.386,  -0.012},
	{  2.08,    0.962,   0.258, -31.2,     0.208,  -0.383,  -6.11,   -0.37,   -0.013},
	{  2.91,    0.974,   0.226, -30.7,     0.23,   -0.375,  -6.64,   -0.453,  -0.024},
	{  2.76,    0.819,   0.344, -27.7,     0.319,  -0.329,  -5.69,   -0.55,    0.05},
	{  2.05,    0.483,   0.362, -28.2,     0.437,  -0.294,  -6,     -0.582,   0.15},
	{  1.5,    0.59,    0.611, -29,      0.68,   -0.23,   -6.2,    -0.595,   0.13},
	{  1.49,    1.21,    0.529, -29.8,     0.1,    -0.21,   -6.4,    -0.637,   0.158},
	{  1.83,   -0.493,   0.298, -38.3,     0.447,  -0.12,   -6.6,    -1.02,    0.24},
	{  1.21,   -1.04,   -2.27,  -35.3,   -0.33,   -0.1,    -6,     -0.84,    0.15 }}

	local s = .2 * alpha
	local k = fix(s)

	if k <= -2 then
		k = -1
	end
	if k >= 9 then
		k = 8
	end
	local da = s - k
	local l = k + fix(1.1 * math.sign(da))
	k = k + 3
	l = l + 3
	local d = {}
	for i=1, 9 do
		d[i-1] = a[k][i] + math.abs(da) * (a[l][i] - a[k][i])
	end
	return d
end

local function f16_model(x, u)
	local xcg = 0.35

	local thtlc, el, ail, rdr = u[0], u[1], u[2], u[3]

	local s = 300
	local b = 30
	local cbar = 11.32
	local rm = 1.57e-3
	local xcgr = .35
	local he = 160.0
	local c1 = -.770
	local c2 = .02755
	local c3 = 1.055e-4
	local c4 = 1.642e-6
	local c5 = .9604
	local c6 = 1.759e-2
	local c7 = 1.792e-5
	local c8 = -.7336
	local c9 = 1.587e-5
	local rtod = 57.29578
	local g = 32.17

	local xd = t_copy(x)
	local vt = x[0]
	local alpha = x[1]*rtod
	local beta = x[2]*rtod
	local phi = x[3]
	local theta = x[4]
	local psi = x[5]
	local p = x[6]
	local q = x[7]
	local r = x[8]
	local alt = x[11]
	local power = x[12]
	
	-- air data computer and engine model
	local amach, qbar = adc(vt, alt)
	local cpow = tgear(thtlc)

	xd[12] = pdot(power, cpow)

	local t = thrust(power, alt, amach)
	local dail = ail / 20
	local drdr = rdr / 30
	
	-- morelli model (polynomial lookup)
	local cxt, cyt, czt, clt, cmt, cnt = morellif16(alpha*math.pi/180, beta*math.pi/180, el*math.pi/180, ail*math.pi/180, rdr*math.pi/180,
		p, q, r, cbar, b, vt, xcg, xcgr)
	
	local tvt = .5 / vt
	local b2v = b * tvt
	local cq = cbar * q * tvt
	
	local d = dampp(alpha)
	local cxt = cxt + cq * d[0]
	local cyt = cyt + b2v * (d[1] * r + d[2] * p)
	local czt = czt + cq * d[3]
	local clt = clt + b2v * (d[4] * r + d[5] * p)
	local cmt = cmt + cq * d[6] + czt * (xcgr-xcg)
	local cnt = cnt + b2v * (d[7] * r + d[8] * p)-cyt * (xcgr-xcg) * cbar/b
	local cbta = math.cos(x[2])
	local u = vt * math.cos(x[1]) * cbta
	local v = vt * math.sin(x[2])
	local w = vt * math.sin(x[1]) * cbta
	local sth = math.sin(theta)
	local cth = math.cos(theta)
	--if math.abs(cth) < 0.01 then
	--	cth = math.sign(cth) * 0.01
	--end
	local sph = math.sin(phi)
	local cph = math.cos(phi)
	local spsi = math.sin(psi)
	local cpsi = math.cos(psi)
	local qs = qbar * s
	local qsb = qs * b
	local rmqs = rm * qs
	local gcth = g * cth
	local qsph = q * sph
	local ay = rmqs * cyt
	local az = rmqs * czt
	
	-- force equations
	local udot = r * v-q * w-g * sth + rm * (qs * cxt + t)
	local vdot = p * w-r * u + gcth * sph + ay
	local wdot = q * u-p * v + gcth * cph + az
	local dum = (u * u + w * w)

	xd[0] = (u * udot + v * vdot + w * wdot)/vt
	xd[1] = (u * wdot-w * udot)/dum
	xd[2] = (vt * vdot-v * xd[0]) * cbta/dum

	-- state kinematics
	xd[3] = p + (sth/cth) * (qsph + r * cph)
	xd[4] = q * cph-r * sph
	xd[5] = (qsph + r * cph)/cth
	-- state moments
	xd[6] = (c2 * p + c1 * r + c4 * he) * q + qsb * (c3 * clt + c4 * cnt)

	xd[7] = (c5 * p-c7 * he) * r + c6 * (r * r-p * p) + qs * cbar * c7 * cmt
	xd[8] = (c8 * p-c2 * r + c9 * he) * q + qsb * (c4 * clt + c9 * cnt)

	-- nav
	local t1 = sph * cpsi
	local t2 = cph * sth
	local t3 = sph * spsi
	local s1 = cth * cpsi
	local s2 = cth * spsi
	local s3 = t1 * sth-cph * spsi
	local s4 = t3 * sth + cph * cpsi
	local s5 = sph * cth
	local s6 = t2 * cpsi + t3
	local s7 = t2 * spsi-t1
	local s8 = cph * cth
	xd[9] = u * s1 + v * s3 + w * s6 -- north speed
	xd[10] = u * s2 + v * s4 + w * s7 -- east speed
	xd[11] = u * sth-v * s5-w * s8 -- vertical speed
	
	local xa = 15.0
	local az = az-xa * xd[7]
	ay = ay+xa*xd[8]
  
	local Nz = (-az / g) - 1
	local Ny = ay / g
	
	return xd, Nz, Ny, az, ay
end

return f16_model
