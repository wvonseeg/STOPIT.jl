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

end

function _printlayer(layer::T) where {T<:AbstractAbsorber}

end

function _printstandardmedia()
    for (i, medium) in enumerate(Absorbers.STANDARDMEDIA)
        println("$i: $medium")
    end
end