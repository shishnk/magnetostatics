var boundariesParameters = BoundaryParameters.ReadJson("input/boundaryParameters.json");
var meshParameters = new MeshParameters("input/meshParameters");
var mesh = new SuperMesh(meshParameters, new LinearMeshBuilder());
var boundaryHandler = new LinearBoundaryHandler(boundariesParameters, meshParameters);
var muB = new MathDependence("mu", "B");
muB.LoadData("input/mu(B)");
var integrator = new Integration(Quadratures.SegmentGaussOrder5());
Spline spline = Spline.CreateBuilder()
    .SetBasis(new HermiteBasis())
    .SetIntegrator(integrator)
    .SetParameters((1E-07, 1E-07))
    .SetPartitions(200)
    .SetPoints(muB.Data!.Select(data => new Point2D(data.Argument, data.Function)).ToArray());
var assembler = new BiMatrixAssembler(new LinearBasis(), integrator, mesh);
SolverFem problem = SolverFem.CreateBuilder()
    .SetMesh(mesh)
    .SetSolverSlae(new CGMCholesky(1000, 1E-15))
    .SetAssembler(assembler)
    .SetDependence(muB)
    .SetNonLinearParameters((1E-10, 100))
    .SetSpline(spline)
    .SetBoundaries(boundaryHandler.Process());

problem.Received += assembler.ReceivePermeability;

problem.Compute();

problem.CalculateAzAtPoint((-.0078, 0.0016));
problem.CalculateAzAtPoint((-.0037, 0.0017));
// problem.CalculateAzAtPoint((-.0013, 0.0019));
// problem.CalculateAzAtPoint((.004, 0.0015));
// problem.CalculateAzAtPoint((.0087, 0.0013));
//
// problem.CalculateBAtPoint((-.0078, 0.0016));
// problem.CalculateBAtPoint((-.0037, 0.0017));
// problem.CalculateBAtPoint((-.0013, 0.0019));
// problem.CalculateBAtPoint((.004, 0.0015));
// problem.CalculateBAtPoint((.0087, 0.0013));