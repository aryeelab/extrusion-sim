import Base.show
using DataFrames
using Random

import Pkg
#Pkg.add("Graphs")
using Graphs

#Pkg.add("Plots")
using Plots

#Pkg.add("ColorSchemes")
using ColorSchemes

using Statistics

mutable struct Cohesin
    processivity::Int64
    alive::Int64
    left_pos::Int64
    right_pos::Int64
    left_blocked::Bool
    right_blocked::Bool
    left_ctcf_captured::Bool
    right_ctcf_captured::Bool
    extruded_beads::Array{Bool}
    
    function Cohesin(processivity, polymer_length, alive=0)
        loading_pos = rand(1:polymer_length)
        new(processivity, alive, loading_pos, loading_pos, false, false, false, false, fill(false, polymer_length))
    end
end
    
function Base.show(io::IO, object::Cohesin) 
        left_blocked_str = string((object.left_blocked) ? " (Cohesin blocked)" : "", 
                                  (object.left_ctcf_captured) ? " (CTCF captured)" : "")
        right_blocked_str = string((object.right_blocked) ? " (Cohesin blocked)" : "", 
                                   (object.right_ctcf_captured) ? " (CTCF captured)" : "") 
        print(io, string("Lived ", object.alive, " of ", object.processivity, " - "), 
                  string("Anchors: ", object.left_pos, left_blocked_str, ", ", object.right_pos, right_blocked_str, " - "), 
                  string("Extruded beads:", findall(object.extruded_beads), "\n")) 
end


mutable struct Sim
    polymer_length::Int64
    bead_size::Int64
    cohesins::Array{Cohesin}
    plus_ctcfs::DataFrame
    minus_ctcfs::DataFrame
    occupied_positions::Array{Bool}
    processivity::Float64
    ctcf_capture_rate::Float64
    ctcf_stabilization_factor::Int64
    extruded_beads::Array{Bool}
    bead_graph::SimpleGraph
    
    function Sim(n, polymer_length_bp=1e6, processivity_bp=150e3, bead_size = 1e3, 
                 ctcf_capture_rate=0.25, ctcf_stabilization_factor=4, plus_ctcfs=missing,
                 minus_ctcfs=missing)
        # Define polymer_length and processity in bead units
        polymer_length = Int(polymer_length_bp/bead_size)
        processivity = processivity_bp / bead_size
        occupied_positions = fill(false, polymer_length)
        alive = 0
        cohesins = []
        bead_graph = path_graph(polymer_length)
        for i in 1:n
            collision = true
            while (collision) 
                global coh = Cohesin(processivity, polymer_length, alive)
                if (!occupied_positions[coh.left_pos])
                    collision = false
                end
            end
            occupied_positions[coh.left_pos] = true
            alive = round(alive + processivity/n)
            push!(cohesins, coh)
            add_edge!(bead_graph, coh.left_pos, copy(coh.left_pos))
        end
        empty_ctcf_df = DataFrame(pos_bp=Float64[], pos=Int64[], occupancy_rate=Float64[], bound_cohesin=Bool[])
        if ismissing(plus_ctcfs)
            plus_ctcfs = deepcopy(empty_ctcf_df)
        end
        if ismissing(minus_ctcfs)
            minus_ctcfs = deepcopy(empty_ctcf_df)
        end
        new(polymer_length, bead_size, cohesins, plus_ctcfs, minus_ctcfs, occupied_positions, 
            processivity, ctcf_capture_rate, ctcf_stabilization_factor, fill(false, polymer_length), bead_graph)            
    end
end

function Base.show(io::IO, object::Sim)
    print(io, string("Cohesin positions: ", findall(object.occupied_positions), "\n"),
              string("Extruded beads: ", sum(object.extruded_beads), " of ", object.polymer_length, "\n"),
              string("Plus CTCF positions (bound fraction): ", object.plus_ctcfs.pos_bp, 
                        "(", object.plus_ctcfs.occupancy_rate, ") ", "\n"),
              string("Minus CTCF positions (bound fraction): ", object.minus_ctcfs.pos_bp, 
                        "(", object.minus_ctcfs.occupancy_rate, ") ", "\n"), "Individual cohesins:")
    for coh in object.cohesins 
        show(coh)
    end
end

function add_ctcf(object::Sim, orientation, pos_bp, occupancy_rate=1) 
    if orientation=="+" 
        push!(object.plus_ctcfs, (pos_bp, ceil(pos_bp/object.bead_size), occupancy_rate, false))
        sort!(object.plus_ctcfs, :pos_bp, rev=true)
    end
    if orientation=="-"
        push!(object.minus_ctcfs, (pos_bp, ceil(pos_bp/object.bead_size), occupancy_rate, false))
        sort!(object.minus_ctcfs, :pos_bp, rev=false)
    end
    object
end

function delete_ctcf(object::Sim, pos_bp)
    idx = indexin(pos_bp, object.plus_ctcfs.pos_bp)[1]
    if idx != nothing
        delete!(object.plus_ctcfs, [idx])
    end
    idx = indexin(pos_bp, object.minus_ctcfs.pos_bp)[1]
    if idx != nothing
        delete!(object.minus_ctcfs, [idx])
    end
    object
end

function advance(object::Sim)
    occupied_positions = object.occupied_positions
    for i in 1:length(object.cohesins)
        coh = object.cohesins[i]
        coh.alive += 1
        if coh.alive >= coh.processivity
            # println("Cohesin fell off!")
            if coh.left_pos+1 != coh.right_pos
                rem_edge!(object.bead_graph, coh.left_pos, coh.right_pos)
            end
            occupied_positions[coh.left_pos] = false
            occupied_positions[coh.right_pos] = false
            # Release CTCF (if any)
            idx = indexin(coh.left_pos, object.plus_ctcfs.pos)[1]
            if idx != nothing
                object.plus_ctcfs.bound_cohesin[idx] = false
            end
            idx = indexin(coh.right_pos, object.minus_ctcfs.pos)[1]
            if idx != nothing
                object.minus_ctcfs.bound_cohesin[idx] = false
            end
            # Generate a new one (making sure it doesn't overlap an existing one)
            collision = true
            while (collision)
                coh = Cohesin(object.processivity, object.polymer_length)
                if (!occupied_positions[coh.left_pos])
                    collision = false
                end
            end
            occupied_positions[coh.left_pos] = true
            add_edge!(object.bead_graph, coh.left_pos, copy(coh.left_pos)) 
        else
            # LEFT
            new_left = coh.left_pos - 1
            if (!coh.left_ctcf_captured)
                # Left side captured by a CTCF?
                for idx in findall((object.plus_ctcfs.pos .== new_left) .& (object.plus_ctcfs.bound_cohesin .== false))
                    if rand(Float64) < object.plus_ctcfs.occupancy_rate[idx] * object.ctcf_capture_rate
                        # println("Captured!")
                        object.plus_ctcfs.bound_cohesin[idx] = true
                        coh.left_ctcf_captured = true
                        coh.processivity = object.processivity * object.ctcf_stabilization_factor
                    end
                end
                if (new_left>=1) 
                    if (occupied_positions[new_left]==false)
                        coh.left_blocked = false
                        if coh.left_pos != coh.right_pos
                            occupied_positions[coh.left_pos] = false
                        end
                        if coh.left_pos+1 != coh.right_pos
                            rem_edge!(object.bead_graph, coh.left_pos, coh.right_pos)
                        end
                        add_edge!(object.bead_graph, new_left, coh.right_pos)
                        coh.left_pos = new_left
                        occupied_positions[new_left] = true
                    end
                else
                    coh.left_blocked = true
                end
            end
            # RIGHT
            new_right = coh.right_pos + 1
            if (!coh.right_ctcf_captured)
                # Right side captured by a CTCF?
                for idx in findall((object.minus_ctcfs.pos .== new_right) .& (object.minus_ctcfs.bound_cohesin .== false))
                    if rand(Float64) < object.minus_ctcfs.occupancy_rate[idx] * object.ctcf_capture_rate
                        # println("Captured!")
                        object.minus_ctcfs.bound_cohesin[idx] = true
                        coh.right_ctcf_captured = true
                        coh.processivity = object.processivity * object.ctcf_stabilization_factor
                    end
                end
                if (new_right<=object.polymer_length)
                    if (occupied_positions[new_right]==false)
                        coh.right_blocked = false
                        if coh.left_pos != coh.right_pos
                            occupied_positions[coh.right_pos] = false
                        end
                        if coh.left_pos+1 != coh.right_pos
                            rem_edge!(object.bead_graph, coh.left_pos, coh.right_pos)
                        end
                        add_edge!(object.bead_graph, coh.left_pos, new_right)
                        coh.right_pos = new_right
                        occupied_positions[new_right] = true
                    end
                else
                    coh.right_blocked = true
                end
            end
            if coh.right_pos - coh.left_pos >= 2
                coh.extruded_beads[((coh.left_pos+1):(coh.right_pos-1))] .= true
            end
        end
        object.cohesins[i] = coh
    end
    object.occupied_positions = occupied_positions
    extruded_beads_sum = fill(false, object.polymer_length)
    for coh in object.cohesins
        extruded_beads_sum = extruded_beads_sum .| coh.extruded_beads
    end
    object.extruded_beads = extruded_beads_sum
    object
end

function effective_dist(object::Sim)
    mat = transpose(dijkstra_shortest_paths(object.bead_graph, 1).dists)
    for i in 2:object.polymer_length
        mat = vcat(mat, transpose(dijkstra_shortest_paths(object.bead_graph, i).dists))
    end
    mat
end

function plot_fiber(object::Sim, xmin=nothing, xmax=nothing, main=nothing)
    x_list = [i for i in 1:object.polymer_length]
    y_list = map(x -> (x ? 2 : 1), object.extruded_beads)
    palette = ColorSchemes.Set1_9
    color_list = [ColorSchemes.grays[1] for i in 1:object.polymer_length]
    shape_list = [:circle for i in 1:object.polymer_length]
    for i in 1:length(object.cohesins) 
        col = palette[i]
        coh = object.cohesins[i]
        color_list[coh.left_pos] = col
        color_list[coh.right_pos] = col
        if coh.left_blocked
            shape_list[coh.left_pos] = :square
        end
        if coh.left_ctcf_captured
            shape_list[coh.left_pos] = :rtriangle
        end
        if coh.right_blocked
            shape_list[coh.right_pos] = :square
        end
        if coh.right_ctcf_captured
            shape_list[coh.right_pos] = :ltriangle
        end
    end
    if isnothing(xmin)
        xmin = (-2)
    end
    if isnothing(xmax)
        xmax = length(object.extruded_beads)+2
    end
    p = scatter(x_list, y_list, color=color_list, markershape=shape_list, label="", markersize=10, markerstrokewidth=0, 
        xlimits=(xmin,xmax), ylimits=(0,5))
    
    for i in 1:length(object.cohesins)
        col = palette[i]
        coh = object.cohesins[i]
        annotate!(mean([coh.left_pos, coh.right_pos]), 2.5+i*0.7, 
            text(string("Cohesin ", i, "\n(", coh.alive, " of ", coh.processivity, ")"), col, :center, 9))
    end
    if ((nrow(object.plus_ctcfs) + nrow(object.minus_ctcfs)) > 0 )
        ctcf_plot_df = [object.plus_ctcfs; object.minus_ctcfs]
        ctcf_plot_df.pch = [[">" for i in 1:nrow(object.plus_ctcfs)]; ["<" for i in 1:nrow(object.minus_ctcfs)]]
        ctcf_plot_df = groupby(ctcf_plot_df, :pos)
        ctcf_plot_df = combine(ctcf_plot_df, :pos, eachindex => :y, :pch)
        for i in 1:nrow(ctcf_plot_df)
            annotate!(ctcf_plot_df.pos[i], 0.25*ctcf_plot_df.y[i], text(ctcf_plot_df.pch[i], :center, 12))
        end
    end
    p
end

default(size = (900, 300))

Random.seed!(1234)
object = Sim(2, 50e3, 10e3, 1000, 0.9, 4)

object = add_ctcf(object, "+", 15000, 1)
object = add_ctcf(object, "+", 15005, 1)
object = add_ctcf(object, "+", 15010, 1)
object = add_ctcf(object, "-", 30000, 1)
object = add_ctcf(object, "-", 30005, 1)
object = add_ctcf(object, "-", 30010, 1)

object

dists = effective_dist(object)
for i in 1:20000
    # plot the first 20 steps in the simulation
    if i<=20
        display(plot_fiber(object))
    end
    advance(object)
    dists = (dists + effective_dist(object))
end
dists = dists/20000

# Average effective distances between points over the duration of the simulation
dists

contacts = (dists.+10).^(-1)

plot(heatmap(dists[10:40,10:40]))

plot(dists[15,:])

plot(heatmap(log.(contacts[10:40,10:40])))


