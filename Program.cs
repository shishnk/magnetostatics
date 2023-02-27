var boundariesParameters = BoundaryParameters.ReadJson("input/boundaryParameters.json");
var meshParameters = new MeshParameters("input/meshParameters");
var mesh = new SuperMesh(meshParameters, new LinearMeshBuilder());
var boundaryHandler = new LinearBoundaryHandler(boundariesParameters, meshParameters);
SolverFem problem = SolverFem.CreateBuilder()
    .SetMesh(mesh)
    // .SetTest(new Test1())
    .SetSolverSlae(new CGMCholesky(1000, 1E-15))
    .SetAssembler(new BiMatrixAssembler(new LinearBasis(), new(Quadratures.SegmentGaussOrder5()), mesh))
    .SetBoundaries(boundaryHandler.Process());

problem.Compute();

// problem.CalculateInPoint(new(-3e-2, 2e-2));