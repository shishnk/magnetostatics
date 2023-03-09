namespace Magnetostatics.src.FEM;

using Spline;

public sealed class FemSolver
{
    public class FemSolverBuilder
    {
        private readonly FemSolver _femSolver = new();

        public FemSolverBuilder SetTest(ITest test)
        {
            _femSolver._test = test;
            return this;
        }

        public FemSolverBuilder SetMesh(IBaseMesh mesh)
        {
            _femSolver._mesh = mesh;
            return this;
        }

        public FemSolverBuilder SetSolverSlae(IterativeSolver iterativeSolver)
        {
            _femSolver._iterativeSolver = iterativeSolver;
            return this;
        }

        public FemSolverBuilder SetBoundaries(IEnumerable<IBoundary> boundaries)
        {
            _femSolver._boundaries = boundaries.DistinctBy(b => b.Node);
            return this;
        }

        public FemSolverBuilder SetAssembler(BaseMatrixAssembler matrixAssembler)
        {
            _femSolver._matrixAssembler = matrixAssembler;
            return this;
        }

        public FemSolverBuilder SetDependence(MathDependence dependence)
        {
            _femSolver._dependence = dependence;
            return this;
        }

        public FemSolverBuilder SetNonLinearParameters((double Residual, int MaxIters) nonLinearParameters)
        {
            _femSolver._nonLinearParameters = nonLinearParameters;
            return this;
        }

        public FemSolverBuilder SetSpline(Spline spline)
        {
            _femSolver._spline = spline;
            return this;
        }

        public static implicit operator FemSolver(FemSolverBuilder builder)
            => builder._femSolver;
    }

    private IBaseMesh _mesh = default!;
    private ITest _test = default!;
    private IterativeSolver _iterativeSolver = default!;
    private IEnumerable<IBoundary> _boundaries = default!;
    private Vector<double> _localVector = default!;
    private Vector<double> _globalVector = default!;
    private BaseMatrixAssembler _matrixAssembler = default!;
    private MathDependence? _dependence;
    private (double Residual, int MaxIters)? _nonLinearParameters;
    private Spline? _spline;
    private bool _isInit;

    public event EventHandler<double>? Received;

    public void Compute()
    {
        Initialize();

        int iter;
        var iters = _nonLinearParameters?.MaxIters ?? 1;

        AssemblySystem();
        AccountingDirichletBoundary();

        var qk = new Vector<double>(_globalVector.Length);

        _isInit = true;

        for (iter = 0; iter < iters; iter++)
        {
            _iterativeSolver.SetMatrix(_matrixAssembler.GlobalMatrix!);
            _iterativeSolver.SetVector(_globalVector);
            _iterativeSolver.Compute();

            qk.Add(_iterativeSolver.Solution!.Value);
            _matrixAssembler.GlobalMatrix!.Clear();
            _globalVector.Fill(0.0);

            AssemblySystem();
            AccountingDirichletBoundary();

            var residual = (_matrixAssembler.GlobalMatrix! * qk - _globalVector).Norm() / _globalVector.Norm();

            Console.WriteLine(residual);

            if (residual < _nonLinearParameters?.Residual) break;
        }

        using var sw = new StreamWriter("output/q.txt");

        foreach (var value in _iterativeSolver.Solution!.Value)
        {
            sw.WriteLine(value);
        }

        Console.WriteLine($"Iterations = {iter}");

        // CalculateError();
    }

    private void OnReceived(double value) => Received?.Invoke(this, value);

    private void Initialize()
    {
        PortraitBuilder.Build(_mesh, out var ig, out var jg);
        _spline?.Compute();
        _matrixAssembler.GlobalMatrix = new(ig.Length - 1, jg.Length)
        {
            Ig = ig,
            Jg = jg
        };

        _globalVector = new(ig.Length - 1);
        _localVector = new(_matrixAssembler.Basis.Size);
    }

    private void AssemblySystem()
    {
        for (int ielem = 0; ielem < _mesh.Elements.Count; ielem++)
        {
            var element = _mesh.Elements[ielem];

            if (!_isInit && element.AreaNumber == 1 && _dependence is not null)
            {
                OnReceived(_dependence.Data![0].Function);
            }

            if (_isInit && element.AreaNumber == 1 && _dependence is not null)
            {
                double mu;
                var localPoints = _mesh.Elements[ielem].Nodes.Select(node => _mesh.Points[node]).ToArray();
                var massCenter = new Point2D(localPoints.Sum(p => p.X) / 4.0, localPoints.Sum(p => p.Y) / 4.0);
                var module = CalculateBAtPoint(massCenter);
                var firstDependenceValue = _dependence!.Data!.First();
                var lastDependenceValue = _dependence.Data!.Last();
                // Console.WriteLine($"|B| = {module}");

                if (module > lastDependenceValue.Argument)
                {
                    mu = 1.0 / (lastDependenceValue.Argument / module * (1.0 / lastDependenceValue.Function - 1.0) +
                                1.0);
                    OnReceived(mu);
                    // Console.WriteLine($"> table value, mu = {mu}");
                }
                else if (module <= firstDependenceValue.Argument)
                {
                    mu = firstDependenceValue.Function;
                    OnReceived(mu);
                    // Console.WriteLine($"< table value, mu = {mu}");
                }
                else
                {
                    OnReceived(_spline!.ValueAtPoint(module));
                    // Console.WriteLine($"inside, mu = {_spline!.ValueAtPoint(module)}");
                }
            }

            _matrixAssembler.BuildLocalMatrices(ielem);
            BuildLocalVector(ielem);

            for (int i = 0; i < _matrixAssembler.Basis.Size; i++)
            {
                _globalVector[element.Nodes[i]] += _localVector[i];

                for (int j = 0; j < _matrixAssembler.Basis.Size; j++)
                {
                    _matrixAssembler.FillGlobalMatrix(element.Nodes[i], element.Nodes[j],
                        _matrixAssembler.StiffnessMatrix[i, j]);
                }
            }
        }
    }

    private void BuildLocalVector(int ielem)
    {
        _localVector.Fill(0.0);
        var jCurrent = _mesh.Areas.First(area => area.Number == _mesh.Elements[ielem].AreaNumber).Current;

        for (int i = 0; i < _matrixAssembler.Basis.Size; i++)
        {
            for (int j = 0; j < _matrixAssembler.Basis.Size; j++)
            {
                // _localVector[i] += _matrixAssembler.MassMatrix[i, j] *
                //                    _test.J(_mesh.Points[_mesh.Elements[ielem].Nodes[j]]);
                _localVector[i] += _matrixAssembler.MassMatrix[i, j] * jCurrent;
            }
        }
    }

    private void AccountingDirichletBoundary()
    {
        // int[] checkBc = new int[_mesh.Points.Count];
        //
        // checkBc.Fill(-1);
        // var boundariesArray = _boundaries.ToArray();
        //
        // // foreach (var b in boundariesArray)
        // // {
        // //     _matrixAssembler.GlobalMatrix.Di[b.Node] = 1E+32;
        // //     _globalVector[b.Node] = 1E+32 * b.Value;
        // // }
        //
        // for (var i = 0; i < boundariesArray.Length; i++)
        // {
        //     boundariesArray[i].Value = _test.Az(_mesh.Points[boundariesArray[i].Node]);
        //     checkBc[boundariesArray[i].Node] = i;
        // }
        //
        // for (int i = 0; i < _mesh.Points.Count; i++)
        // {
        //     int index;
        //     if (checkBc[i] != -1)
        //     {
        //         _matrixAssembler.GlobalMatrix!.Di[i] = 1.0;
        //         _globalVector[i] = boundariesArray[checkBc[i]].Value;
        //
        //         for (int k = _matrixAssembler.GlobalMatrix.Ig[i]; k < _matrixAssembler.GlobalMatrix.Ig[i + 1]; k++)
        //         {
        //             index = _matrixAssembler.GlobalMatrix.Jg[k];
        //
        //             if (checkBc[index] == -1)
        //             {
        //                 _globalVector[index] -= _matrixAssembler.GlobalMatrix.Gg[k] * _globalVector[i];
        //             }
        //
        //             _matrixAssembler.GlobalMatrix.Gg[k] = 0.0;
        //         }
        //     }
        //     else
        //     {
        //         for (int k = _matrixAssembler.GlobalMatrix!.Ig[i]; k < _matrixAssembler.GlobalMatrix.Ig[i + 1]; k++)
        //         {
        //             index = _matrixAssembler.GlobalMatrix.Jg[k];
        //
        //             if (checkBc[index] == -1) continue;
        //             _globalVector[i] -= _matrixAssembler.GlobalMatrix.Gg[k] * _globalVector[index];
        //             _matrixAssembler.GlobalMatrix.Gg[k] = 0.0;
        //         }
        //     }
        // }

        var boundariesArray = _boundaries.ToArray();

        foreach (var boundary in boundariesArray)
        {
            _matrixAssembler.GlobalMatrix!.Di[boundary.Node] = 1E+32;
            _globalVector[boundary.Node] = 0.0;
        }
    }

    private void CalculateError()
    {
        var error = new double[_mesh.Points.Count];

        for (int i = 0; i < error.Length; i++)
        {
            error[i] = Math.Abs(_iterativeSolver.Solution!.Value[i] - _test.Az(_mesh.Points[i]));
        }

        Array.ForEach(error, Console.WriteLine);

        var sum = error.Sum(t => t * t);

        sum = Math.Sqrt(sum / _mesh.Points.Count);

        Console.WriteLine($"rms = {sum}");

        // using var sw = new StreamWriter("output/3.csv");
        //
        // for (int i = 0; i < error.Length; i++)
        // {
        //     if (i == 0)
        //     {
        //         sw.WriteLine("$i$, $u_i^*$, $u_i$, $|u^* - u|$, Погрешность");
        //         sw.WriteLine(
        //             $"{i}, {_test.U(_mesh.Points[i])}, {_iterativeSolver.Solution!.Value[i]}, {error[i]}, {sum}");
        //         continue;
        //     }
        //
        //     sw.WriteLine($"{i}, {_test.U(_mesh.Points[i])}, {_iterativeSolver.Solution!.Value[i]}, {error[i]},");
        // }
    }

    public double CalculateAzAtPoint(Point2D point)
    {
        var res = 0.0;

        var ielem = FindElementNumber(point);

        var element = _mesh.Elements[ielem];
        var bPoint = _mesh.Points[element.Nodes[0]];
        var ePoint = _mesh.Points[element.Nodes[^1]];

        double hx = ePoint.X - bPoint.X;
        double hy = ePoint.Y - bPoint.Y;

        var ksi = (point.X - bPoint.X) / hx;
        var eta = (point.Y - bPoint.Y) / hy;

        var templatePoint = new Point2D(ksi, eta);

        for (int i = 0; i < _matrixAssembler.Basis.Size; i++)
        {
            res += _iterativeSolver.Solution!.Value[_mesh.Elements[ielem].Nodes[i]] *
                   _matrixAssembler.Basis.GetPsi(i, templatePoint);
        }

        Console.WriteLine($"Az at {point} = {res}");
        return res;
    }

    public double CalculateBAtPoint(Point2D point)
    {
        var ielem = FindElementNumber(point);
        OnReceived(1.0 / PhysicsConstants.VacuumPermeability);
        _matrixAssembler.BuildLocalMatrices(ielem);

        var sqrModule = 0.0;

        for (int i = 0; i < _matrixAssembler.Basis.Size; i++)
        {
            for (int j = 0; j < _matrixAssembler.Basis.Size; j++)
            {
                sqrModule += _matrixAssembler.StiffnessMatrix[i, j] *
                             _iterativeSolver.Solution!.Value[_mesh.Elements[ielem].Nodes[i]] *
                             _iterativeSolver.Solution.Value[_mesh.Elements[ielem].Nodes[j]];
            }
        }

        var elementArea = (_mesh.Points[_mesh.Elements[ielem].Nodes[1]].X -
                           _mesh.Points[_mesh.Elements[ielem].Nodes[0]].X) *
                          (_mesh.Points[_mesh.Elements[ielem].Nodes[2]].Y -
                           _mesh.Points[_mesh.Elements[ielem].Nodes[0]].Y);

        sqrModule /= elementArea;

        var module = Math.Sqrt(sqrModule);
        //
        // var dx = _iterativeSolver.Solution!.Value[_mesh.Elements[ielem].Nodes[1]] -
        //          _iterativeSolver.Solution!.Value[_mesh.Elements[ielem].Nodes[0]];
        // var dy = _iterativeSolver.Solution!.Value[_mesh.Elements[ielem].Nodes[2]] -
        //          _iterativeSolver.Solution!.Value[_mesh.Elements[ielem].Nodes[0]];
        //
        // var element = _mesh.Elements[ielem];
        // var bPoint = _mesh.Points[element.Nodes[0]];
        // var ePoint = _mesh.Points[element.Nodes[^1]];
        //
        // double hx = ePoint.X - bPoint.X;
        // double hy = ePoint.Y - bPoint.Y;
        //
        // var rotor = (dy / hy, -dx / hx);

        // Console.WriteLine($"|B| at ({point.X}; {point.Y}) = {module}");
        // Console.WriteLine(
        // $"|B| at {point} = {rotor.Item1} {rotor.Item2} = {Math.Sqrt(rotor.Item1 * rotor.Item1 + rotor.Item2 * rotor.Item2)}");
        return module;
        // return Math.Sqrt(rotor.Item1 * rotor.Item1 + rotor.Item2 * rotor.Item2);

        // var ielem = FindElementNumber(point);
        //
        // var element = _mesh.Elements[ielem];
        // var bPoint = _mesh.Points[element.Nodes[0]];
        // var ePoint = _mesh.Points[element.Nodes[^1]];
        //
        // double hx = ePoint.X - bPoint.X;
        // double hy = ePoint.Y - bPoint.Y;
        //
        // var ksi = (point.X - bPoint.X) / hx;
        // var eta = (point.Y - bPoint.Y) / hy;
        //
        // var templatePoint = (ksi, eta);
        //
        // double dx = 0.0;
        // double dy = 0.0;
        //
        // for (int i = 0; i < _matrixAssembler.Basis.Size; i++)
        // {
        //     dx += _iterativeSolver.Solution!.Value[element.Nodes[i]] *
        //           _matrixAssembler.Basis.GetDPsi(i, 0, templatePoint);
        //     dy += _iterativeSolver.Solution!.Value[element.Nodes[i]] *
        //           _matrixAssembler.Basis.GetDPsi(i, 1, templatePoint);
        // }
        //
        // dx = -dx / hx;
        // dy /= hy;
        //
        // return Math.Sqrt(dx * dx + dy * dy); // calculate with B_x and B_y components of rotor A
    }

    private int FindElementNumber(Point2D point)
    {
        foreach (var (element, idx) in _mesh.Elements.Select((element, idx) => (element, idx)))
        {
            if (point.X >= _mesh.Points[element.Nodes[0]].X && point.X <= _mesh.Points[element.Nodes[1]].X &&
                point.Y >= _mesh.Points[element.Nodes[0]].Y && point.Y <= _mesh.Points[element.Nodes[2]].Y)
            {
                return idx;
            }
        }

        throw new("Not supported exception!");
    }

    public static FemSolverBuilder CreateBuilder() => new();
}