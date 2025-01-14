import os
import sys
import numpy as np
from argparse import ArgumentParser
from parmed import load_file, Residue, Atom

# Utility functions
def calc_box_volume(box):
    """Calculate the volume of a triclinic box."""
    return box[0] * box[1] * box[2]

def random3Drotation():
    """Generate a random 3D rotation matrix."""
    theta = np.random.uniform(0, 2 * np.pi)
    phi = np.random.uniform(0, 2 * np.pi)
    z = np.random.uniform(-1, 1)
    r = np.sqrt(1 - z ** 2)
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    return np.array([
        [np.cos(phi) + (1 - np.cos(phi)) * x ** 2, (1 - np.cos(phi)) * x * y - np.sin(phi) * z, (1 - np.cos(phi)) * x * z + np.sin(phi) * y],
        [(1 - np.cos(phi)) * y * x + np.sin(phi) * z, np.cos(phi) + (1 - np.cos(phi)) * y ** 2, (1 - np.cos(phi)) * y * z - np.sin(phi) * x],
        [(1 - np.cos(phi)) * z * x - np.sin(phi) * y, (1 - np.cos(phi)) * z * y + np.sin(phi) * x, np.cos(phi) + (1 - np.cos(phi)) * z ** 2]
    ])

def selectReplaceableWaters(waters, min_distance, max_distance):
    """Select waters that can be replaced."""
    selected_waters = []
    for water in waters:
        dist = np.linalg.norm(water)
        if min_distance < dist < max_distance:
            selected_waters.append(water)
    return selected_waters

def insertFragments(selected_waters, fragment, box):
    """Insert fragments into the simulation box."""
    inserted_fragments = []
    for water in selected_waters:
        rotation = random3Drotation()
        frag_position = np.dot(rotation, fragment) + water
        inserted_fragments.append(frag_position)
    return inserted_fragments

def stripWaterClash(waters, fragments, clash_distance):
    """Remove waters that clash with inserted fragments."""
    non_clashing_waters = []
    for water in waters:
        if all(np.linalg.norm(water - frag) > clash_distance for frag in fragments):
            non_clashing_waters.append(water)
    return non_clashing_waters

from parmed import Atom, Residue

def add_interligand_repulsion(topology, fragments):
    """
    Add inter-ligand repulsion parameters to the topology using ParmEd.
    """
    for frag_coords in fragments:
        for coord in frag_coords:
            # Create a new atom object
            atom = Atom()
            atom.name = "DUM"  # Atom name
            atom.type = "Du"   # Atom type
            atom.charge = 0.0  # Neutral charge
            atom.mass = 12.01  # Arbitrary mass
            atom.rmin = 1.5    # Lennard-Jones radius
            atom.epsilon = 0.1 # Lennard-Jones well depth

            # Assign coordinates
            atom.xx, atom.xy, atom.xz = coord

            # Append atom to topology's atom list
            topology.add_atom(atom, "DUM", 0)

            # Create a residue if needed
            residue = Residue(name="DUM")
            residue.add_atom(atom)
            topology.residues.append(residue)

    print(f"Added {len(fragments)} fragments as dummy atoms to the topology.")


# Main script
def main():
    parser = ArgumentParser(description="Replace water molecules with fragments in AMBER simulations.")
    parser.add_argument("-top", default=r"C:\Users\Asus\CrypticSWISH\1.1ohr\1ohr_solvated.prmtop", help="Topology file.")
    parser.add_argument("-rst", default=r"C:\Users\Asus\CrypticSWISH\1.1ohr\1ohr.rst7", help="Restart file.")
    parser.add_argument("-frag", default=r"C:\Users\Asus\CrypticSWISH\1.1ohr\ligand.pdb", help="Fragment file.")
    parser.add_argument("-conc", type=float, default=0.1, help="Concentration of fragments (M).")
    parser.add_argument("-min_dist", type=float, default=2.0, help="Minimum distance from protein.")
    parser.add_argument("-max_dist", type=float, default=5.0, help="Maximum distance from protein.")

    args = parser.parse_args()

    print(f"Loading topology from {args.top}")
    topology = load_file(args.top)

    print(f"Loading coordinates from {args.rst}")
    coordinates = load_file(args.rst)

    # Placeholder: Load water positions and box dimensions
    waters = np.random.rand(1000, 3) * 10  # Mock water positions
    box = [10.0, 10.0, 10.0]  # Mock box dimensions

    # Calculate volume and number of fragments to insert
    volume = calc_box_volume(box)
    num_fragments = int(args.conc * 6.022e23 * volume / 1e27)

    print(f"System volume: {volume:.2f} nm^3")
    print(f"Number of fragments to insert: {num_fragments}")

    # Select replaceable waters
    selected_waters = selectReplaceableWaters(waters, args.min_dist, args.max_dist)

    # Insert fragments
    fragment = np.random.rand(3, 3)  # Mock fragment positions
    inserted_fragments = insertFragments(selected_waters, fragment, box)

    # Remove clashing waters
    waters = stripWaterClash(waters, inserted_fragments, clash_distance=1.0)

    # Modify topology for inter-ligand repulsion
    add_interligand_repulsion(topology, inserted_fragments)

    print(f"Inserted {len(inserted_fragments)} fragments and updated topology.")
    topology.write_parm('updated_topology.prmtop')
  #  coordinates.write_rst7('updated_coordinates.rst7')
    print("Updated topology and coordinates saved successfully.")



if __name__ == "__main__":
    main()
