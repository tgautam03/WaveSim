function wave_sim_1d(x_details::Dict{String, Int64}, t_details::Dict{String, Int64}, src_details::Dict{String, Any}, c::Float64, op::Int64)
    # Space Discretization
	x_min = x_details["x_min"]
	x_max = x_details["x_max"]
	nx = x_details["nx"]
	dx = (x_max-x_min)/(nx-1) 
	x = range(start=x_min, stop=x_max, length=nx)

	# Time Discretization
	t_min = t_details["t_min"]
	t_max = t_details["t_max"]
	nt = t_details["nt"]
	dt = (t_max-t_min)/(nt-1)
	t = range(start=t_min, stop=t_max, length=nt)
    
    # Source Function
	f0 = src_details["f0"] # Dominant Frequency
	t0 = src_details["t0"] # Source Time Shift
	isrc = Int(floor(src_details["isrc"]/dx)+1) # Source Location
	
	src = src_details["src"] # Source Function

    ########################################################
    ############ Finite Difference Simulation ##############
    ########################################################
    # Pressure Field at p(x,t)
	p = zeros(nx)
	# Pressure Field at p(x,t-dt)
	p_prev = zeros(nx)
	# Pressure Field at p(x,t+dt)
	p_next = zeros(nx)
	# Solution at each step
	p_sols = zeros(nt,nx)

    # Looping over time
    for it in range(start=2, stop=nt)
        if op == 3
            # Looping over Space
            for ix in range(start=2, stop=nx-1)
                # Evaluating 2nd derivative wrt x
                d2p_dx2 = p[ix+1] - 2*p[ix] + p[ix-1]
                # Updating Solution
                if ix == isrc
                    p_next[ix] = (c*dt/dx)^2 * d2p_dx2 + 2*p[ix] - p_prev[ix] + dt^2 * src[it]
                else
                    p_next[ix] = (c*dt/dx)^2 * d2p_dx2 + 2*p[ix] - p_prev[ix]
                end
            end

        elseif op == 5
            # Looping over Space
            for ix in range(start=3, stop=nx-2)
                # Evaluating 2nd derivative wrt x
                d2p_dx2 = -1/12 * p[ix+2] + 4/3  * p[ix+1] - 5/2 * p[ix] +4/3  * p[ix - 1] - 1/12 * p[ix - 2]
                # Updating Solution
                if ix == isrc
                    p_next[ix] = (c*dt/dx)^2 * d2p_dx2 + 2*p[ix] - p_prev[ix] + dt^2 * src[it]
                else
                    p_next[ix] = (c*dt/dx)^2 * d2p_dx2 + 2*p[ix] - p_prev[ix]
                end
            end
        end

        

        # Current Sol becomes Previous Sol
        p_prev[:] = p[:]
        # Next Sol becomes Current Sol
        p[:] = p_next[:]
    
        # Storing solutions at each time step
        p_sols[it,:] = p
    end

    return p_sols
end