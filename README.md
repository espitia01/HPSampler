# Sampling the conformal space of an HP protein using Monte Carlo Methods

<img src = "loop.png">

## Monte Carlo Loop Steps

### Proposed Move
- **Description**: The initial step where a new conformation of the protein is proposed. This is usually done by slightly modifying the current conformation.

### Random Position
- **Description**: A position is chosen at random in the sequence for potential movement or change.

### Initialization
- **Description**: The initial setup of the sequence and its conformation before the Monte Carlo simulation begins.

### Record Data
- **Description**: At this step, relevant data about the current state of the sequence is recorded, which might include energy, position, and other attributes.

### Overlap Check
- **Description**: Before accepting the proposed move, it's necessary to ensure that the new conformation doesn't result in any overlaps (i.e., two amino acids occupying the same position in space).

### Energy Evaluation
- **Description**: The energy of the new conformation is evaluated. This is often based on the interactions between hydrophobic (H) and polar (P) amino acids.

### Update Optimal Configuration
- **Description**: If the new conformation is deemed better (usually having lower energy), then it's accepted as the new optimal configuration.

### Metropolis Criteria
- **Description**: A probabilistic criterion used to decide whether to accept or reject the proposed move. Even if the new conformation has higher energy, it might still be accepted with a certain probability to allow for exploration of the conformational space.

## Hashing Mechanism in HP Protein Folding Simulation

The simulation uses a hashing mechanism to efficiently track seen protein configurations. This is done to handle the isosymmetry property of the protein configurations. In the 2D plane, a protein's conformation might be identical to another when subjected to operations like translation, rotation, or reflection. These symmetrical conformations are equivalent in terms of their structural properties and energies.

### Handling Isosymmetry:

1. **Translation**: The `translate_to_origin` function ensures that a conformation is translated such that its starting point is at the origin.

2. **Rotation**: The protein configuration can be rotated by 0, 90, 180, or 270 degrees. The `rotate` function handles this transformation using a basic rotation matrix.

3. **Reflection**: A protein configuration can also be reflected over the y-axis, which is managed by the `reflect` function.

After these transformations, the `canonical_form` function determines the canonical (or standard) form of a given protein conformation by finding the 'smallest' configuration lexicographically among all possible transformed forms. 

### Hashing:

The `hash_conformation` function produces a unique hash for a protein's canonical form. This is done by:
- Transforming the protein conformation to its canonical form.
- Flattening this form into a 1D tuple.
- Generating a hash value for this tuple.

Throughout the Monte Carlo loop, each time a new protein conformation is proposed, its hash is calculated. If this hash matches one in the set of previously seen hashes (`seen_hashes`), the configuration is considered already explored. The simulation can then skip redundant evaluations, increasing efficiency.

This approach ensures that the simulation doesn't repeatedly evaluate energetically identical configurations, thus speeding up the search for the minimum energy state while reducing computational overhead.

