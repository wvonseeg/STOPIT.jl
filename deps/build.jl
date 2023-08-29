using Artifacts
import Pkg

artifactTOML = find_artifacts_toml(@__FILE__)

# Query the `Artifacts.toml` file for an existing hash
AME2020hash = artifact_hash("AME2020", artifactTOML)

# If the name was not bound, or the hash it was bound to does not exist, create it!
if AME2020hash === nothing || !artifact_exists(AME2020hash)
    # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
    AME2020hash = Pkg.Artifacts.create_artifact() do artifact_dir
        # We create the artifact by simply downloading files into the artifact directory
        download("https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt", joinpath(artifact_dir, "ame2020.txt"))
    end

    # Now bind that hash within our `Artifacts.toml`
    Pkg.Artifacts.bind_artifact!(artifactTOML, "AME2020", AME2020hash; force=true)
end

Pkg.ensure_artifact_installed("AME2020", artifactTOML)