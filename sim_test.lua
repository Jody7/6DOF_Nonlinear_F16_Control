local sim_f16 = require("sim_f16.lua")
local ctrl_f16 = require("ctrl_f16")
local RK4 = require("RK4")

local function t_copy(original)
	local copy = {}
	for key, value in pairs(original) do
		copy[key] = value
	end
	return copy
end

local frame_counter = 0
local function visualize_plane(x)
--visualize with your game engine/platform of choice
end

local function print_state(_state)
	local state = t_copy(_state)
	local out_string = ""
	for i=0, 15 do
		if (i == 3) or (i == 4) or (i == 5) then
			state[i] = math.round(100*math.deg(state[i]))/100 .. " DEG"
		end
		print("state[" .. i .. "]: " .. tostring(state[i]))
	end
	print("----")
end

local x_state = {}
for i=0, 15 do
	x_state[i] = 0.001
end
local u_input = {}
for i=0, 3 do
	u_input[i] = 0
end

x_state[0] = 400
x_state[5] = math.rad(0)
x_state[11] = 4000

u_input[0] = 0
u_input[1] = 0
u_input[2] = 0
u_input[3] = 0
_G.saved_x_state = t_copy(x_state)

local left_yaw_amt = 0
local right_yaw_amt = 0
local throttle_delta_amt = 0
local throttle_amt = 0

--TODO: for those referencing with repo, you simply need to integrate the game engine/control interface of your choice
-- to input into u_input vector

dt = (0.001)
local g_t = 0
while true do
	throttle_delta_amt = math.min(1, math.max(-1, throttle_delta_amt))
	throttle_amt = math.min(math.max(-0.1, throttle_amt + (throttle_delta_amt * dt)), 1.0)
	if frame_counter % 30 == 0 then
		print("dt:", dt)
		print_state(x_state)
	end
	local display_u_deg = nil
	local display_Nz = nil
	local display_ps = nil
	local display_Ny_r = nil
  
	u_input[3] = throttle_amt
	frame_counter = frame_counter + 1
	
	local rk4_subdivision = 8
	for i=1, rk4_subdivision do	
		local temp_t, temp_x = RK4.rk4(x_state, function(t, x_hat)
			--local x_d_hat = sim_f16(x_hat, u_input)
			local x_d_hat_ctrl_state, u_deg, Nz, ps, Ny_r = ctrl_f16(x_hat, u_input, sim_f16)
			
			display_u_deg = u_deg
			display_Nz = Nz
			display_ps = ps
			display_Ny_r = Ny_r
			
			--return x_d_hat
			return x_d_hat_ctrl_state
		end, 0, dt / rk4_subdivision)
		x_state = temp_x
		
		x_state[0] = math.min(10000, math.max(0, x_state[0]))
		x_state[1] = math.min(2, math.max(-2, x_state[1]))
		x_state[2] = math.min(2, math.max(-2, x_state[2]))
		x_state[3] = math.min(99, math.max(-99, x_state[3]))
		x_state[4] = math.min(99, math.max(-99, x_state[4]))
		x_state[5] = math.min(99, math.max(-99, x_state[5]))
		x_state[6] = math.min(3, math.max(-3, x_state[6]))
		x_state[7] = math.min(3, math.max(-3, x_state[7]))
		x_state[8] = math.min(3, math.max(-3, x_state[8]))
		--x_state[9] = math.min(10000, math.max(-10000, x_state[9]))
		--x_state[10] = math.min(10000, math.max(-10000, x_state[10]))
		x_state[11] = math.min(100000, math.max(-1000, x_state[11]))
		x_state[12] = math.min(100, math.max(0, x_state[12]))
	end
	
	dt = wait()
end

