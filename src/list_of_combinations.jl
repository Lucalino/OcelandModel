# This function takes two arrays and outputs a new array which's elements
# are arrays with each combination of the elements of the two input lists.
# It can be used to sample points in the parameter space.
using IterTools

function listofpoints(arrays)
    arr = Array{Float64}[]
    for p in product(arrays...)
        push!(arr,[y for y in p])
    end
    return arr
end

# function listofpoints(list1::Array{Float64},list2::Array{Float64})
#     arr = Array{Float64}[]
#     for p in product(list1,list2)
#         push!(arr,[y for y in p])
#     end
#     return arr
# end    
