## Voronoi Game 🎨
Julia implementation of a two-player game, played on a [Voronoi tessellation](https://en.wikipedia.org/wiki/Voronoi_diagram#Mathematics), using [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation). Part of the class "Computerorientierte Mathematik II" @ TU Berlin (Summer Term 2025).

## Usage
The folder containing "src" must be called "Voronoi" with a capital V!
Setup the project by running this command inside the Voronoi folder:
```
julia> using Pkg; Pkg.activate("."); Pkg.instantiate()
```
Then the project can be started with this command:
```
julia> using Pkg; Pkg.activate("."); using Voronoi; main()
```