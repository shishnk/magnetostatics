namespace Magnetostatics.src.Spline;

public class Spline
{
    public class SplineBuilder
    {
        private readonly Spline _spline = new();

        public SplineBuilder SetParameters((double Alpha, double Beta) parameters)
        {
            _spline._parameters = parameters;
            return this;
        }

        public SplineBuilder SetPoints(Point2D[] points)
        {
            _spline._points = points;
            return this;
        }

        public SplineBuilder SetPartitions(int partitions)
        {
            _spline._partitions = partitions;
            return this;
        }

        public SplineBuilder SetBasis(IBasis1D basis)
        {
            _spline._basis = basis;
            return this;
        }

        public SplineBuilder SetIntegrator(Integration integrator)
        {
            _spline._integrator = integrator;
            return this;
        }

        public static implicit operator Spline(SplineBuilder builder)
            => builder._spline;
    }

    private Interval[] _elements = default!;
    private Point2D[] _points = default!;
    private SevenDiagonalMatrix _globalMatrix = default!;
    private Matrix _localMatrix = default!;
    private FEM.Vector<double> _vector = default!;
    private List<Point2D> _result = default!;
    private IBasis1D _basis = default!;
    private Integration _integrator = default!;
    private (double Alpha, double Beta) _parameters;
    private int _partitions;

    public IReadOnlyList<Point2D> InitialPoints => _points;
    public IReadOnlyList<Point2D> Result => _result;

    private void Init()
    {
        _globalMatrix = new(_elements.Length * 2 + 2);
        _localMatrix = new(_basis.Size);
        _vector = new(_globalMatrix.Size);
        _result = new();
    }

    private void BuildMesh()
    {
        _elements = new Interval[_partitions];
        _points = _points.OrderBy(p => p.X).ToArray();

        if (_partitions == 1)
        {
            _elements[0] = new(_points[0].X, _points[^1].X);
            return;
        }

        double step = (_points.MaxBy(p => p.X)!.X - _points.MinBy(p => p.X)!.X) / _partitions;
        _elements[0] = new(_points[0].X, _points[0].X + step);

        for (int i = 1; i < _elements.Length; i++)
        {
            _elements[i] = new(_elements[i - 1].RightBorder, _elements[i - 1].RightBorder + step);
        }
    }

    public void Compute()
    {
        BuildMesh();
        Init();
        AssemblyMatrix();
        var gaussSeidel = new GaussSeidel(1000, 1E-15, 1.23);
        _vector = gaussSeidel.Compute(_globalMatrix, _vector);
        // ValuesAtPoints();
    }

    public double ValueAtPoint(double point)
    {
        int ielem = -1;
        double result = 0.0;

        for (int i = 0; i < _elements.Length; i++)
        {
            if (!_elements[i].IsContain(point)) continue;
            ielem = i;
            break;
        }

        if (ielem == -1) throw new("Not supported exception!");

        double x = (point - _elements[ielem].LeftBorder) / _elements[ielem].Length;

        for (int i = 0; i < _basis.Size; i++)
        {
            result += _vector[2 * ielem + i] * _basis.GetPsi(i, x, _elements[ielem].Length);
        }

        return result;
    }

    private void AssemblyMatrix()
    {
        int[] checker = new int[_points.Length];
        checker.Fill(1);

        for (int ielem = 0; ielem < _elements.Length; ielem++)
        {
            for (int ipoint = 0; ipoint < _points.Length; ipoint++)
            {
                if (!_elements[ielem].IsContain(_points[ipoint].X) || checker[ipoint] != 1) continue;

                checker[ipoint] = -1;
                double x = (_points[ipoint].X - _elements[ielem].LeftBorder) / _elements[ielem].Length;

                for (int i = 0; i < _basis.Size; i++)
                {
                    _vector[2 * ielem + i] += _points[ipoint].Y * _basis.GetPsi(i, x, _elements[ielem].Length);

                    for (int j = 0; j < _basis.Size; j++)
                    {
                        double Function1(double point, double h)
                        {
                            var dFi1 = _basis.GetDPsi(i, point, h);
                            var dFi2 = _basis.GetDPsi(j, point, h);

                            return dFi1 * dFi2;
                        }

                        double Function2(double point, double h)
                        {
                            var ddFi1 = _basis.GetDdPsi(i, point, h);
                            var ddFi2 = _basis.GetDdPsi(j, point, h);

                            return ddFi1 * ddFi2;
                        }

                        _localMatrix[i, j] += _basis.GetPsi(i, x, _elements[ielem].Length) *
                                              _basis.GetPsi(j, x, _elements[ielem].Length) +
                                              _parameters.Alpha * _integrator.Gauss1D(Function1, _elements[ielem]) +
                                              _parameters.Beta * _integrator.Gauss1D(Function2, _elements[ielem]);
                    }
                }
            }

            for (int i = 0; i < _localMatrix.Size; i++)
            {
                _globalMatrix[0, 2 * ielem + i] += _localMatrix[i, i];
            }

            for (int i = 0; i < _localMatrix.Size - 1; i++)
            {
                _globalMatrix[1, 2 * ielem + i] += _localMatrix[i, i + 1];
                _globalMatrix[4, 2 * ielem + i] += _localMatrix[i, i + 1];
            }

            for (int i = 0; i < _localMatrix.Size - 2; i++)
            {
                _globalMatrix[2, 2 * ielem + i] = _localMatrix[i, i + 2];
                _globalMatrix[5, 2 * ielem + i] = _localMatrix[i, i + 2];
            }

            for (int i = 0; i < _localMatrix.Size - 3; i++)
            {
                _globalMatrix[3, 2 * ielem + i] = _localMatrix[i, i + 3];
                _globalMatrix[6, 2 * ielem + i] = _localMatrix[i, i + 3];
            }

            _localMatrix.Clear();
        }
    }

    private void ValuesAtPoints()
    {
        double sum = 0.0;

        for (int ielem = 0; ielem < _elements.Length; ielem++)
        {
            Point2D changedPoint = new(_elements[ielem].LeftBorder, 0.0);

            do
            {
                double x = (changedPoint.X - _elements[ielem].LeftBorder) / _elements[ielem].Length;

                for (int i = 0; i < _basis.Size; i++)
                {
                    sum += _vector[2 * ielem + i] * _basis.GetPsi(i, x, _elements[ielem].Length);
                }

                _result.Add(changedPoint with { Y = sum });

                changedPoint += (0.0005, 0.0);
                sum = 0.0;
            } while (_elements[ielem].IsContain(changedPoint.X));
        }
    }

    public static SplineBuilder CreateBuilder()
        => new();
}