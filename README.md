# AMS 562 Final Project - Spherical-Triangular-Meshes #

In this project as the final project in AMS 562(Scientific Programming in C++ in Stony Brook University), we implement geometric calculation on spherical
triangle meshes.

A mesh (unstructured) can be represented by a set of coordinates (geometric data)
and a connectivity table (topological data) that forms the connection relation
between triangles and coordinates.

There're three algorithms that declared in `src/manip.hpp`:
 - Determine node to triangle adjacency
 - Compute averaged outward unit normal vectors
 - Analyze the approximation errors

More instructions are listed in final-proj.pdf.

# Manual to compile
To compile
```console
make
```

To do testing
```console
make test
```

To clean directory
```console
make clean
```

To debug with `valgrind`, type

```console
make debug
```

To get the result of test after compile
```console
./main
```
