function _printmenu()
    println("\n1 Define stoppee")
    println("2 Define absorber")
    println("3 Edit absorber")
    println("4 Run with current parameters")
    println("5 Find thickness of absorber to stop the stopee")
    println("6 Print status of data")
    println("7 Stop\n")
end

function _printhelp()

end

function _printstatus()
    println("Stopee => ")
    display(stopee)
    for (i, layer) in enumerate(sandwich.layers)
        println("Layer $i")
        display(layer)
    end
end

function _printstandardmedia()
    println("\nStandard Media")
    for (i, medium) in enumerate(Absorbers.STANDARDMEDIA)
        println("$i: $medium")
    end
end

function _printresults(energyloss::Vector{typeof(1.0u"MeV")})
    print("Stopee => ")
    display(stopee)
    for (i, layer) in enumerate(sandwich.layers)
        println("Layer $i")
        display(layer)
        println("Energy lost = $(energyloss[i])")
    end
    remaining = stopee.energy - sum(energyloss)
    println("\nEnergy remaining = $(remaining)")
end