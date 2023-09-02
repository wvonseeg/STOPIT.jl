module Interactive

using ..STOPIT, ..STOPIT.Absorbers, ..STOPIT.Particles
using NaturallyUnitful
import ReplMaker: initrepl

export initsession

global stopee::Particle
global sandwich::Sandwich = Sandwich()


include("printing.jl")

function _parseinput(input::String)
    # General commands
    if input == "?"
        _printhelp()
        return
    elseif input == "m"
        _printmenu()
        return
    end

    # Menu commands
    if input == "1"
        println("Define stopee")
    elseif input == "2"
        println("Define absorber")
    elseif input == "3"
        println("Edit absorber")
    elseif input == "4"
        println("Run with current parameters")
    elseif input == "5"
        println("Find stopping thickness for stopee")
    elseif input == "6"
        println("Print status")
    elseif input == "7"
        println("Quit")
    else
        println("Command not recognized!")
    end

    _printmenu()
end

function initsession()
    initrepl(_parseinput, prompt_text="STOPIT> ", prompt_color=:blue,
        start_key=')', mode_name="STOPIT Mode")

    println("Type \"?\" for help or \"m\" to see the menu")
end
end