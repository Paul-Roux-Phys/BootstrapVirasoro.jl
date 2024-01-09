#==================

CFTdata.jl contains a module CFTdata that provides types representing
central charges and fields in 2D CFTs with Virasoro symmetry.

==================#


"""
Provides types representing central charges and fields in CFT.
"""
module CFTData

using Match;

export CentralCharge,Field

#print complex numbers in latex format
function Base.show(io::IO,::MIME"text/latex",z::Complex)
    print(real(z)," + ",imag(z),"i")
end

#=

Central charges

=#

#=  The central charge can be given directly or as a function of b,β,B. No matter what
    parameter is passed we compute B from it, then b,β,c. 
    We ensure that the relation between b and beta never changes (sqrt branch). 
=#

"""Get B from given parameter"""
Bfrom(parameter,value) = @match parameter begin
        "c" => (value-13+sqrt(complex((value-1)*(value-25))))/12
        "b" => value^2
        "β" => -value^2
        "B" => value
    end

"""Get asked parameter from B"""
Bto(parameter,value) = @match parameter begin
        "c" => 13+6*value+6/value
        "b" => sqrt(complex(value))
        "β" => im*sqrt(complex(value))
        "B" => value
    end

"""
    CentralCharge{T}
Object representing the central charge. 
Contains the values of the 4 parameters representing it.
"""
struct CentralCharge{T}

    #= T is the type of the parameters; either Complex{Float64} or Complex{BigFloat} 
       for arbitrary precision. =#
    values::Dict{String,T}
    
end

"""
    CentralCharge(parameter,value)

Constructor function for the CentralCharge type.

Given one of the four parameters `"c"`,`"b"`,`"β"`,`"B"` and its value, 
creates an object CentralCharge{T} where T is the type of `value`.

#Example 
```julia-repl
julia> setprecision(BigFloat,20,base=10)
julia> charge = CentralCharge("β",sqrt(big(2)))
Central charge :
B = -2.0 + 0.0im
c = -2.0 + 0.0im
b = 0.0 + 1.414213562373095048804im
β = -1.414213562373095048804 + 0.0im
```
"""
function CentralCharge(parameter="c",value=1)
    # Constructor
    T=typeof(AbstractFloat(real(value)))
    B=Bfrom(parameter,value)
    dict=Dict(key => Bto(key,B) for key in ("c","b","β","B"))
    CentralCharge{complex(T)}(dict)
end

"""Display an object of type CentralCharge"""
function Base.show(io::IO,charge::CentralCharge{T}) where {T}
    println("Central charge :")
    for (key,value) in charge.values
        println(io, "$key = $value")
    end
end

"""Display the value of the central charge in LaTeX format"""
function Base.show(io::IO,::MIME"text/latex",charge::CentralCharge{T},parameter) where {T}
    #Latex copy-pastable output
    if parameter=="β"
        print("\\beta = ")
    else
        print(parameter," = ")
    end
    show(io,MIME("text/latex"),charge[parameter])
end

"""Overload of [] to access values in charge"""
Base.getindex(charge::CentralCharge,key) = charge.values[key];

#=

Fields

=#

"""Get p from any given parameter"""
p_from(parameter,value,charge::CentralCharge) = @match parameter begin
        "Δ" => sqrt(complex(value - (charge["c"]-1)/24))
        "δ" => sqrt(complex(value))
        "P" => -im*value
        "p" => value
    end

"""Get all parameters from p"""
p_to(parameter,value,charge::CentralCharge) = @match parameter begin
        "Δ" => value^2 + (charge["c"]-1)/24
        "δ" => value^2
        "P" => im*value
        "p" => value
    end

"""
    Field{T}
Object representing a conformal field. 
Contains the values of the 4 parameters `"Δ"`,`"δ"`,`"P"`,`"p"` for its conformal dimension,
and flags saying whether the field is in the Kac table, degenerate, logarithmic or diagonal.
"""
struct Field{T}

    values::Dict{String,Vector{T}}
    isKac::Bool
    r::Rational
    s::Rational
    isdegenerate::Bool
    islogarithmic::Bool
    isdiagonal::Bool

end

"""
    Field(charge, parameter, leftvalue, rightvalue; kwargs...)

Constructor function for the Field type.

Given a charge `charge`, one of the four parameters `"Δ"`,`"δ"`,`"P"`,`"p"` and two values, 
create an object Field{T} (where T is the type of the values in `charge`) that represents a 
field of left and right dimensions given by leftvalue and rightvalue in the chosen 
parametrisation.

# keyword arguments: 

* `Kac::Bool`: if set to true, the field can be constructed from the values of its r and s 
indices,
* `r::Rational`,`s::Rational`: used in conjunction to `Kac=true`, must be given rational 
values,
* `logarithmic::Bool`: set to True if the field is logarithmic,
* `degenerate::Bool`: set to True if the field is degenerate,
* `diagonal::Bool`: set to True to get a diagonal field ; only the leftvalue needs to be 
given.

#Examples
```julia-repl
julia> charge = CentralCharge("b",big(0.5));
julia> field = Field(charge,Kac=true,r=0,s=1)
Non-diagonal field with Kac indices r = 0//1, s = 1//1 and (left,right) dimensions:
Δ = ( 2.5625 + 0.0im , 2.5625 + 0.0im )
P = ( -0.0 - 1.0im , 0.0 + 1.0im )
δ = ( 1.0 - 0.0im , 1.0 + 0.0im )
p = ( -1.0 + 0.0im , 1.0 + 0.0im )
```
```julia-repl
julia> charge = CentralCharge("β",1.5+im);
julia> Field(charge,"δ",2,3)
Non-diagonal field with (left,right) dimensions:
Δ = ( 2.1579142011834325 - 0.6789940828402367im , 3.1579142011834316 - 0.6789940828402367im )
P = ( 0.0 + 1.4142135623730951im , 0.0 + 1.7320508075688772im )
δ = ( 2.0000000000000004 + 0.0im , 2.9999999999999996 + 0.0im )
p = ( 1.4142135623730951 + 0.0im , 1.7320508075688772 + 0.0im )
```
```julia-repl
julia> charge = CentralCharge();
julia> Field(charge,"δ",1,diagonal=true)
Diagonal field of dimension:
Δ = 1.0 + 0.0im
P = 0.0 + 1.0im
δ = 1.0 + 0.0im
p = 1.0 + 0.0im
```
"""
function Field(
    charge::CentralCharge = CentralCharge("c",1),parameter="Δ",leftvalue=0,rightvalue=0;
    Kac=false,r=0,s=0,logarithmic=false,degenerate=false,diagonal=false
)
    #Constructor
    T=typeof(charge.values["c"])   #dimensions have the same type as central charges
    if degenerate
        Kac=true
    end
    if Kac
        pleft = 1/2*(charge["b"]*r - 1/charge["b"]*s)
        pright = 1/2*(charge["b"]*r + 1/charge["b"]*s)
    else
        pleft,pright = p_from.(parameter,[leftvalue,rightvalue],Ref(charge))
    end
    if diagonal
        pright=pleft
    end
    values = Dict(key => p_to.(key,[pleft,pright],Ref(charge)) for key in ("Δ","δ","P","p"))
    Field{complex(T)}(values,Kac,r,s,degenerate,logarithmic,diagonal)
end

"""Compute the spin Δleft - Δright of a field."""
function spin(field::Field)
    #Computes the spin Δ-Δbar
    if field.isdiagonal
        return 0
    else
        return field["Δ"][1]-field["Δ"][2]
    end
end

"""Display field"""
function Base.show(io::IO,field::Field)
    #Print fields
    if field.isdiagonal
        println("Diagonal field of dimension:")
        for (key,value) in field.values
            println(io, key," = ", value[1])
        end
    else
        print("Non-diagonal field")
        if field.isKac
            print(" with Kac indices r = ",field.r,", s = ",field.s, " and ")
        else
            print(" with ")
        end
        println("(left,right) dimensions:")
        for (key,value) in field.values
                println(io, key," = ( ", value[1]," , ",value[2]," )")
        end
    end
end

"""Display dimension of field in latex format"""
function Base.show(io::IO,::MIME"text/latex",field::Field{T},parameter) where {T}
    if field.isdiagonal
        if parameter=="Δ"
            print("\\Delta = ")
        elseif parameter=="δ"
            print("\\delta = ")
        else
            print(parameter," = ")
        end
        show(io,MIME("text/latex"),field[parameter][1])
    else
        if parameter=="Δ"
            print("(\\Delta,\\bar\\Delta) = ")
        elseif parameter=="δ"
            print("(\\delta,\\bar\\delta) = ")
        else
            print("(",parameter,"\\bar",parameter,") = ")
        end
        print("("),show(io,MIME("text/latex"),field[parameter][1]);print(",");
        show(io,MIME("text/latex"),field[parameter][2]);print(")")
    end
end

"""Overload []"""
Base.getindex(field::Field,key) = field.values[key];

end