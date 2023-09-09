function _editabsorber()
    len = length(sandwich.layers)
    println("Which layer [1-$len]?")
    index = parse(Int, readline())
    if index < 1 || index > len
        println("Layer $index does not exist")
        return
    end
    _definelayer(index)
end

function _defineabsorber()
    println("How many layers in the absorber?")
    numlayers = parse(Int, readline())
    for i in 1:numlayers
        _definelayer(i)
    end
end

function _definelayer(layernum::Integer)
    println("##### Layer $layernum #####\n")
    println("1: Solid")
    println("2: Gas")
    println("3: Standard medium")
    med = parse(Int, readline())
    if med == 1
        _solidmedium(layernum)
    elseif med == 2
        _gaseousmedium(layernum)
    elseif med == 3
        _standardmedium(layernum)
    end
end

function _solidmedium(layernum::Integer)
    println("Enter density [g/cm^3] and thickness [mg/cm^2]")
    dens = 1.0u"g/cm^3" * parse(Float64, readline())
    thick = 1.0u"mg/cm^2" * parse(Float64, readline())
    println("How many elements in layer $layernum?")
    numel = parse(Int, readline())
    A = Vector{<:Integer}(undef, numel)
    Z = Vector{<:Integer}(undef, numel)
    num = Vector{<:Integer}(undef, numel)
    for j in 1:numel
        println("For Layer $layernum, Element $j enter A, Z, and number of atoms per molecule")
        A[j] = parse(Int, readline())
        Z[j] = parse(Int, readline())
        num[j] = parse(Int, readline())
    end
    usrlayer = SolidAbsorber(A, Z, num, thick, dens)
    setlayer!(sandwich, layernum, usrlayer)
end

function _gaseousmedium(layernum::Integer)
    println("Enter pressure [Torr] and depth [cm]")
    pres = 1.0u"Torr" * parse(Float64, readline())
    dep = 1.0u"cm" * parse(Float64, readline())
    println("How many elements in layer $layernum?")
    numel = parse(Int, readline())
    A = Vector{<:Integer}(undef, numel)
    Z = Vector{<:Integer}(undef, numel)
    num = Vector{<:Integer}(undef, numel)
    for j in 1:numel
        println("For Layer $layernum, Element $j enter A, Z, and number of atoms per molecule")
        A[j] = parse(Int, readline())
        Z[j] = parse(Int, readline())
        num[j] = parse(Int, readline())
    end
    usrlayer = GasAbsorber(A, Z, num, pres, dep)
    setlayer!(sandwich, layernum, usrlayer)
end

function _standardmedium(layernum::Integer)
    _printstandardmedia()
    medium = parse(Int, readline())
    pressure = 0.0u"Torr"
    if medium in (1, 4, 5, 8, 9, 11, 12)
        # Gases
        println("Pressure [Torr]?")
        pressure = 1.0u"Torr" * parse(Float64, readline())
    end
    println("Depth [cm]?")
    depth = 1.0u"cm" * parse(Float64, readline())
    usrlayer = getstandardmedium(Absorbers.STANDARDMEDIA[medium]; pressure=pressure, depth=depth)
    setlayer!(sandwich, layernum, usrlayer)
end

function _definestopee()
    println("What is the stopee A?")
    A = parse(Int, readline())
    println("What is the stopee Z?")
    Z = parse(Int, readline())
    println("What is the particle energy [MeV]?")
    E = parse(Float64, readline())
    println("\nStopee Defined:")
    try
        global stopee = Particle(A, Z, E * 1.0u"MeV")
    catch err
        if isa(err, MassNotFoundException)
            println("The requested particle ($A$(STOPIT.atomicsymbols[Z+1])) does not exist.")
            println("Stopee not defined!")
        else
            rethrow()
        end
    end
end

function _run()
    lost = desorb(stopee, sandwich)
    _printresults(lost)
end

function _findthickness()
    println("Find stopping thickness for stopee")
end