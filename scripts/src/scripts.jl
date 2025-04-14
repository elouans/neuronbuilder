module scripts

include("minimal_hh_single.jl")
include("full_hh_single.jl")

greet() = print("Hello World!")


function something()
    println("BBBB")
    full_HH()
    println("BBefefBB")
end

function check_files()
    println("Current directory: ", pwd())
    println("Does full_hh_single.jl exist? ", isfile("full_hh_single.jl"))
    println("Does minimal_hh_single.jl exist? ", isfile("minimal_hh_single.jl"))
end


export check_file, full_HH, minimal_HH
end # module scripts
