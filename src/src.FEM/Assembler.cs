namespace Magnetostatics.src.FEM;

public abstract class BaseMatrixAssembler
{
    protected readonly IBaseMesh _mesh;
    protected readonly Integration _integrator;
    protected Matrix[]? _baseStiffnessMatrix;
    protected Matrix? _baseMassMatrix;
    protected double? _mu;

    public SparseMatrix? GlobalMatrix { get; set; } // need initialize with portrait builder 
    public Matrix StiffnessMatrix { get; }
    public Matrix MassMatrix { get; }
    public IBasis2D Basis { get; }

    protected BaseMatrixAssembler(IBasis2D basis, Integration integrator, IBaseMesh mesh)
    {
        Basis = basis;
        _integrator = integrator;
        _mesh = mesh;
        StiffnessMatrix = new(basis.Size);
        MassMatrix = new(basis.Size);
    }

    public abstract void BuildLocalMatrices(int ielem);

    public void FillGlobalMatrix(int i, int j, double value)
    {
        if (GlobalMatrix is null)
        {
            throw new("Initialize the global matrix (use portrait builder)!");
        }

        if (i == j)
        {
            GlobalMatrix.Di[i] += value;
            return;
        }

        if (i <= j) return;
        for (int ind = GlobalMatrix.Ig[i]; ind < GlobalMatrix.Ig[i + 1]; ind++)
        {
            if (GlobalMatrix.Jg[ind] != j) continue;
            GlobalMatrix.Gg[ind] += value;
            return;
        }
    }
}

public class BiMatrixAssembler : BaseMatrixAssembler
{
    public BiMatrixAssembler(IBasis2D basis, Integration integrator, IBaseMesh mesh) : base(basis, integrator, mesh)
    {
    }

    public override void BuildLocalMatrices(int ielem)
    {
        var element = _mesh.Elements[ielem];
        var bPoint = _mesh.Points[element.Nodes[0]];
        var ePoint = _mesh.Points[element.Nodes[^1]];

        double hx = ePoint.X - bPoint.X;
        double hy = ePoint.Y - bPoint.Y;

        if (_baseStiffnessMatrix is null)
        {
            _baseStiffnessMatrix = new Matrix[] { new(Basis.Size), new(Basis.Size) };
            _baseMassMatrix = new(Basis.Size);
            var templateElement = new Rectangle(new(0.0, 0.0), new(1.0, 1.0));

            for (int i = 0; i < Basis.Size; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    Func<Point2D, double> function;

                    for (int k = 0; k < 2; k++)
                    {
                        var ik = i;
                        var jk = j;
                        var k1 = k;
                        function = p =>
                        {
                            var dFi1 = Basis.GetDPsi(ik, k1, p);
                            var dFi2 = Basis.GetDPsi(jk, k1, p);

                            return dFi1 * dFi2;
                        };

                        _baseStiffnessMatrix[k][i, j] = _baseStiffnessMatrix[k][j, i] =
                            _integrator.Gauss2D(function, templateElement);
                    }

                    var i1 = i;
                    var j1 = j;
                    function = p =>
                    {
                        var fi1 = Basis.GetPsi(i1, p);
                        var fi2 = Basis.GetPsi(j1, p);

                        return fi1 * fi2;
                    };
                    _baseMassMatrix[i, j] = _baseMassMatrix[j, i] = _integrator.Gauss2D(function, templateElement);
                }
            }
        }

        var mu = _mu ?? _mesh.Areas.First(area => area.Number == _mesh.Elements[ielem].AreaNumber).Permeability;

        for (int i = 0; i < Basis.Size; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                StiffnessMatrix[i, j] = StiffnessMatrix[j, i] =
                    1.0 / mu * (hy / hx * _baseStiffnessMatrix[0][i, j] +
                                hx / hy * _baseStiffnessMatrix[1][i, j]);
            }
        }

        _mu = null;

        for (int i = 0; i < Basis.Size; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                MassMatrix[i, j] = MassMatrix[j, i] = hx * hy * _baseMassMatrix![i, j];
            }
        }
    }

    public void ReceivePermeability(object? sender, double value) => _mu = PhysicsConstants.VacuumPermeability * value;
}