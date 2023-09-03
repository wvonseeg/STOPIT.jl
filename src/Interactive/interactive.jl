module Interactive

using ..STOPIT, ..STOPIT.Absorbers, ..STOPIT.Particles
using NaturallyUnitful
import ReplMaker: initrepl

export initsession

global stopee::Particle
global sandwich::Sandwich = Sandwich()

include("printing.jl")
include("functions.jl")

function _parseinput(input::String)
    if input == "?"
        _printhelp()
    elseif input == "m"
        _printmenu()
    elseif input == "1"
        _definestopee()
    elseif input == "2"
        _defineabsorber()
    elseif input == "3"
        _editabsorber()
    elseif input == "4"
        _run()
    elseif input == "5"
        _findthickness()
    elseif input == "6"
        _printstatus()
    elseif input == "7"
        println("Quit")
    else
        println("Command not recognized!")
    end
end

function initsession()
    initrepl(_parseinput, prompt_text="STOPIT> ", prompt_color=:blue,
        start_key=')', mode_name="STOPIT Mode")

    println("Type \"?\" for help or \"m\" to see the menu")

    _printmenu()
end
end