"""Script to create a heater HDF5 fixture for the FlowBoiling IO test."""

import h5py
from types import SimpleNamespace

def createHeaterfile(heaterInfo):
    """
    Create heater file
    """
    filename = h5py.File(heaterInfo["name"] + "_hdf5_htr_0001", "w")

    nsites = heaterInfo["numSites"]
    xsite = [0.0] * nsites
    ysite = [1.0e-13] * nsites
    zsite = [0.0] * nsites
    radii = [0.2] * nsites

    if nsites > 1:
        xspan = heaterInfo["xmax"] - heaterInfo["xmin"]
        zspan = heaterInfo["zmax"] - heaterInfo["zmin"]
        for i in range(nsites):
            frac = i / (nsites - 1)
            xsite[i] = heaterInfo["xmin"] + frac * xspan
            if heaterInfo["dim"] == 1:
                zsite[i] = 0.0
            elif heaterInfo["dim"] == 2:
                zsite[i] = heaterInfo["zmin"] + frac * zspan
            else:
                raise ValueError("Error in Heater.dim")

    filename.create_dataset(
        "heater/xMin", data=heaterInfo["xmin"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/xMax", data=heaterInfo["xmax"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/zMin", data=heaterInfo["zmin"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/zMax", data=heaterInfo["zmax"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/yMin", data=heaterInfo["ymin"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/yMax", data=heaterInfo["ymax"], shape=(1), dtype="float32"
    )

    filename.create_dataset(
        "heater/wallTemp", data=heaterInfo["wallTemp"], shape=(1), dtype="float32"
    )

    filename.create_dataset(
        "heater/advAngle", data=heaterInfo["advAngle"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/rcdAngle", data=heaterInfo["rcdAngle"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/velContact", data=heaterInfo["velContact"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/nucWaitTime", data=heaterInfo["nucWaitTime"], shape=(1), dtype="float32"
    )

    filename.create_dataset(
        "site/num", data=heaterInfo["numSites"], shape=(1), dtype="int32"
    )
    filename.create_dataset(
        "site/x", data=xsite, shape=(heaterInfo["numSites"]), dtype="float32"
    )
    filename.create_dataset(
        "site/y", data=ysite, shape=(heaterInfo["numSites"]), dtype="float32"
    )
    filename.create_dataset(
        "site/z", data=zsite, shape=(heaterInfo["numSites"]), dtype="float32"
    )

    filename.create_dataset(
        "init/radii", data=radii, shape=(heaterInfo["numSites"]), dtype="float32"
    )

    filename.close()

    print(f"Wrote heater information to file {filename}")


if __name__ == "__main__":

    Heater = SimpleNamespace()
  
    # Heater name and number of heaters
    Heater.name = "flow_boiling"
    Heater.dim = 1
    Heater.xmin = -24.0
    Heater.xmax = 24.0
    Heater.zmin = -2.5
    Heater.zmax = 2.5
    Heater.ymin = 0.0
    Heater.ymax = 1e-13

    # Wall temperature
    Heater.wallTemp = 1.0

    # Contact angle information
    Heater.advAngle = 90.0
    Heater.rcdAngle = 45.0
    Heater.velContact = 0.2
    Heater.nucWaitTime = 0.4
    Heater.nucSeedRadius = 0.1

    # Number of nucleation sites
    Heater.numSites = 200

    createHeaterfile(vars(Heater))
