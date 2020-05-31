module TSPImages

using LightGraphs
using Images
using Random
using TravelingSalesmanHeuristics: solve_tsp
using Gadfly
using LinearAlgebra: norm

#=
Note: throughout this package I index pixels in images according to the indexing of normal
Julia matrices. That is, index (1,1) is the top left corner, (1,2) is directly to the right
of the top left corner, and (2,1) is directly below the top left corner. So pixel indices
are of the form (y,x) where increasing y corresponds to moving down. Note that this
indexing convention doesn't match that of many image editing programs.
=#

"Load an image to grayscale. If `invert`, invert the grayscale image. Return a tuple
(image_grayscale, image)."
function load_gray_image(path; invert = false)
    img = load(path)
    img_gray = Gray.(img)
    if invert
        img_gray .= oneunit(eltype(img_gray)) .- img_gray
    end
    return img_gray, img
end

"Sample `n` points according to a grayscale image. If the image has height and width
`(h,w)`, points returned will be in (0,h] x (0,w]. The probability of sampling a given
point is proportional to the darkness of the corresponding pixel in the input image.

Return a 2 by n matrix with one column per sampled point"
function sample_from_image(img::Matrix, n::Int)
    h, w = size(img)
    points = zeros(2,n)
    for pt_ind in 1:n
        # we sample with (1 - rand()) to get values in (0, 1]
        while true
            pt_y = h * (1 - rand())
            pt_x = w * (1 - rand())
            pt_pixel_val = img[ceil(Int, pt_y), ceil(Int, pt_x)]
            if rand() >= pt_pixel_val # accept this point, darker => more likely to accept
                points[:, pt_ind] .= [pt_y, pt_x]
                break
            end
        end
    end
    return points
end

"Partition a matrix of points into components. Two points in separate components are at
least `dist_cutoff` pixels from each other.

Internally, this function constructs a simple graph over the passed points. This graph
has an edge between any two points with distance less than `dist_cutoff`, and we then
return the connected components of this graph.

Return a Vector{Vector{Int}} where each `i` in `1:size(pts, 2)` appears exactly once."
function get_components(pts, dist_cutoff)
    graph = LightGraphs.euclidean_graph(pts; cutoff = dist_cutoff)[1]
    return LightGraphs.connected_components(graph)
end

"Order a matrix of points according to an approximate TSP solution. This will increase
the number of columns in the `pts` matrix by one as the path start/end will be duplicated.

Optional kwargs are `quality_factor` which is passed to `solve_tsp` and `dist_cutoff`
which specifies the maximum allowable edge distance. A value of nothing corresponds to
no maximum."
function order_points(pts; quality_factor = 50, dist_cutoff = nothing)
    n = size(pts, 2)
    distmat = [norm(pts[:,i] - pts[:,j]) for i in 1:n, j in 1:n]
    if !(isnothing(dist_cutoff))
        distmat[distmat .> dist_cutoff] .*= 100
    end
    path, _ = solve_tsp(distmat; quality_factor = quality_factor)
    return pts[:, path]
end

"Create a Gadfly plot layer from ordered points. `pts` should be a 2 by n matrix where
columns correspond to points and the order of columns defines a TSP path."
function points_to_layer(pts, colors; styling = (), frac = 1.0)
    n = size(pts, 2)
    nplot = 1 + round(Int, frac * (n-1))
    return layer(
        x = pts[2, 1:(nplot-1)],
        xend = pts[2, 2:nplot],
        y = -pts[1, 1:(nplot-1)],
        yend = -pts[1, 2:nplot],
        color = colors,
        Geom.segment,
        style(;styling...)
    )
end

"Take a collection of layers and build a plot which combines the layers."
function build_plot(layers, height, width)
    p = plot(layers[1], Coord.cartesian(xmin = 0, xmax = width, ymin = -height, ymax = 0))
    for i in 2:length(layers)
        push!(p, layers[i])
    end
    return p
end

"Take a plot, set the background to white, and remove labels/grids."
function _whiten_plot!(p)
    push!(p, style(background_color = colorant"white"))
    push!(p, Guide.xlabel(nothing))
    push!(p, Guide.ylabel(nothing))
    push!(p, Guide.xticks(;ticks = nothing))
    push!(p, Guide.yticks(;ticks = nothing))
    return nothing
end

"Load an image, sample points, solve TSP subproblems, and plot!

Positional arguments:
- `path`: path to image file
- `n::Int`: number of points to sample from the image
- `component_dist::Real`: distance in pixels between image components. If a group of sapmled
  points is separated from all other sampled points by at least `component_dist` pixels,
  then that group of points will form its own TSP path.

Keyword arguments:
- `layer_styling`
- `quality_factor::Real`: quality factor keyword passed to `solve_tsp`. See
  `TravelingSalesmanHeuristics.solve_tsp` for more information. Defaults to 50.
- `invert_image::Bool`: Whether to invert the grayscale input image so that light areas are
  sampled more than dark areas.
- `path_dist_cutoff::Union{Real, Nothing}`: maximum distance in pixels of an edge in any
  TSP path, or `nothing` for no limit.
- `verbose::Bool`: whether to print timing information during execution.

Return a NamedTuple with keys :plot, :height, and :width corresponding to
the plot object and dimensions of the input image."
function image_to_tsp_plot(
    path, n, component_dist;
    layer_styling = (:default_color => colorant"black",),
    quality_factor = 50,
    invert_image = false,
    path_dist_cutoff = nothing,
    single_color = nothing,
    verbose = true
)
    _time_load = @elapsed img_gray, img = load_gray_image(path; invert = invert_image)
    verbose && @show _time_load
    height, width = size(img)
    _time_sample = @elapsed pts = sample_from_image(img_gray, n)
    verbose && @show _time_sample
    _time_components = @elapsed components = get_components(pts, component_dist)
    verbose && @show _time_components
    _time_tsp = @elapsed layers = map(components) do comp
        pts_ordered = order_points(pts[:, comp];
                                   quality_factor = quality_factor,
                                   dist_cutoff = path_dist_cutoff)
        segment_midpoints = (pts_ordered[:,1:(end-1)] .+ pts_ordered[:,2:end]) / 2
        if isnothing(single_color)
            colors = [img[ceil(Int, col[1]), ceil(Int, col[2])] for col in eachcol(segment_midpoints)]
        else
            colors = [single_color]
        end
        return points_to_layer(pts_ordered, colors; styling = layer_styling)
    end
    verbose && @show _time_tsp
    p = build_plot(layers, height, width)
    _whiten_plot!(p)
    return (plot = p, width = width, height = height)
end

function save_plot(outfile, plot, width, height)
    draw(SVG(outfile, 1px * width, 1px * height), plot)
    return nothing
end

end # module
