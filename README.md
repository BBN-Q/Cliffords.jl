# Cliffords

This library allows for efficient calculation of Clifford circuits by tracking the evolution of X and Z generators (the so-called **tableau** representation). No special effort has been made to strictly minimize the number of bits needed to store each Clifford. Rather, the goal was clarity. One unique feature compared to other such utilities is that we also efficiently track the inverse operations. This is useful to, e.g., compute 'undo' gates in randomized benchmarking sequences.
