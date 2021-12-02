local function t_copy(original)
	local copy = {}
	for key, value in pairs(original) do
		copy[key] = value
	end
	return copy
end

local xqeuil = {502.0, 0.0389, 0.0, 0.0, 0.0389, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1000.0, 9.0567}
local uequil = {0.1395, -0.7496, 0.0, 0.0}
local neg_K_lqr = 
{{156.880151,  31.037008,  38.729833,  -0,        -0,        -0, -0,        -0,      },
{ -0,        -0,        -0,       -37.84483,   25.40956,    6.82876, 332.88343,   17.15997 },
{ -0,        -0,        -0,        23.91233,   -5.69968,   21.63431, -64.4949,    88.36203 }}

local function get_u_deg(u_ref, x_state)
	local x_delta = t_copy(x_state)
	for i=1, #xqeuil do
		x_delta[i-1] = x_delta[i-1] - xqeuil[i]
	end
	local x_ctrl = {}
	x_ctrl[0] = x_state[1]
	x_ctrl[1] = x_state[7]
	x_ctrl[2] = x_state[13]
	x_ctrl[3] = x_state[2]
	x_ctrl[4] = x_state[6]
	x_ctrl[5] = x_state[8]
	x_ctrl[6] = x_state[14]
	x_ctrl[7] = x_state[15]
	
	local u_deg = {}
	for i=0, 3 do
		u_deg[i] = 0
	end
	
	for row_i=1, 3 do
		local dot_sum = 0
		for col_i=1, 8 do
			dot_sum = dot_sum + (x_ctrl[col_i-1]*neg_K_lqr[row_i][col_i])
		end
		u_deg[row_i] = dot_sum
	end
	
	u_deg[0] = u_ref[3]
	for i=0, 3 do
		u_deg[i] = u_deg[i] + uequil[i+1]
	end
	
	-- crappy anti integral
	if math.abs(u_ref[1]) < 0.05 then
		x_ctrl[6] = 0
	end
	
	-- clamp control inputs
	u_deg[0] = math.max(math.min(u_deg[0], 1), 0)
	u_deg[1] = math.max(math.min(u_deg[1], 25), -25)
	u_deg[2] = math.max(math.min(u_deg[2], 21), -21)
	u_deg[3] = math.max(math.min(u_deg[3], 30), -30)
	
	return x_ctrl, u_deg
end

local function get_integrator_derivs(u_ref, Nz, ps, Ny_r)
	local drv = {Nz - u_ref[0], ps - u_ref[1], Ny_r - u_ref[2]}
	return drv
end

local f_counter = 0
local function ctrl_f16(x_state, u_ref, sim_f16)
	--f_counter = f_counter + 1
	local x_0_13 = {}
	for i=0, 12 do
		x_0_13[i] = x_state[i]
	end
	
	local x_ctrl, u_deg = get_u_deg(u_ref, x_state)
	
	--if f_counter % 30 == 0 then
	--	for i=0, 3 do
	--		print(i, u_deg[i])
	--	end
	--	print("===")
	--end
	
	local xd_model, Nz, Ny, _, _ = sim_f16(x_0_13, u_deg)
	local ps = x_ctrl[4] * math.cos(x_ctrl[0]) + x_ctrl[5] * math.sin(x_ctrl[0])
	local Ny_r = Ny + x_ctrl[5]
	
	local xd = {}
	local integrator_derivs = get_integrator_derivs(u_ref, Nz, ps, Ny_r)
	for i=13, 15 do
		xd[i] = integrator_derivs[i-13 + 1]
	end
	for i=0, 12 do
		xd[i] = xd_model[i]
	end

	return xd, u_deg, Nz, ps, Ny_r
end

return ctrl_f16
