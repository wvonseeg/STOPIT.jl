testsandwich = Sandwich()
silayer = getstandardmedium("Si"; depth=0.05u"cm")
proton = Particle(1, 1, 1.0, 10u"MeV")
setlayer!(testsandwich, silayer)

desorb!(proton, testsandwich)