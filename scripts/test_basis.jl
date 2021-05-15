using FiniteElements

h = 1/3
grid = Grid(h)

support = Set([Element(Set([grid.points[6], grid.points[1], grid.points[5]])),
                Element(Set([grid.points[6], grid.points[1], grid.points[2]])),
                Element(Set([grid.points[6], grid.points[2], grid.points[7]])),
                Element(Set([grid.points[6], grid.points[7], grid.points[11]])),
                Element(Set([grid.points[6], grid.points[7], grid.points[10]])),
                Element(Set([grid.points[6], grid.points[5], grid.points[10]]))])

basis1 = BasisFunction(support)
