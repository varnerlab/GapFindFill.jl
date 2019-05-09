#=
Names w/o apostrophe
function 1: read name lists, e.g. reaction names, compound names
  --> return an array of strings
function 2: read bounds, e.g. upper/lower bounds
  --> return an array of tuples (compound::String, bound::Float64)
function 3: read S matrix
  --> return an array of tuples
      (compound::String, reaction::String, coefficient::Float64)
Q: can be combined together? Yup
=#

function read_GAMS_format_file(filePath::String)
  returnArr = Array{Any, 1}()
  open(filePath, "r") do file
    lines = readlines(file)[2:(end-1)]
    if (length(split(lines[1])) <= 1)
      # return an array of strings
      read_names_list(lines, returnArr)
    else
      tokens = split(split(lines[1])[1], ".")
      if length(tokens) == 1
        read_bounds_list(lines, returnArr)
      else
        read_S_matrix_list(lines, returnArr)
      end
    end
  end

  return returnArr
end

function read_names_list(rawData::Array, newData::Array)
  for (id, val) in enumerate(rawData)
    val = split(val)
    if length(val) == 0
      # println("empty line found")
    else
      push!(newData, val[1][2: (end-1)])
    end
  end
end

function read_bounds_list(rawData::Array, newData::Array)
  for (id, val) in enumerate(rawData)
    tokens = split(val)
    if length(tokens) < 1
      # println("weird line found")
    else
      push!(newData, (tokens[1][2:end-1], parse(Float64, tokens[2])))
    end
  end
end

function read_S_matrix_list(rawData::Array, newData::Array)
  for (id, val) in enumerate(rawData)
    tk1 = split(val)
    if length(tk1) == 0
      # println("empty line found")
    else
      tk2 = split(tk1[1], ".")
      push!(newData, (tk2[1][2:end-1], tk2[2][2:end-1], parse(Float64, tk1[2])))
    end
  end
end
