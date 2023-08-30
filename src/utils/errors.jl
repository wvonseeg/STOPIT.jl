export ZeroEnergyException, ParticleStoppedException, MassNotFoundException

struct ZeroEnergyException <: Exception end

struct ParticleStoppedException <: Exception end

struct MassNotFoundException <: Exception end