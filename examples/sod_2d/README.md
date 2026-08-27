# Rotated 2D Sod validation

This experiment embeds the Sod shock tube in a Cartesian 2D grid, first with
the discontinuity normal to x and then normal to y. The two solutions should
match after transposing the grid and exchanging velocity components.

Run the multi-resolution comparison from the repository root:

```bash
python benchmarks/convergence/rotated_sod_2d.py
```

The script also evolves a uniform periodic flow and reports the largest state
change and global conservation errors before producing the rotated-Sod figure.
