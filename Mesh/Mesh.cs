namespace Magnetostatics.Mesh;

public interface IBaseMesh
{
    IReadOnlyList<Point2D> Points { get; }
    IReadOnlyList<FiniteElement> Elements { get; }
}

public abstract class MeshBuilder
{
    protected abstract int ElementSize { get; }

    public abstract (IEnumerable<Point2D>, FiniteElement[]) Build(MeshParameters meshParameters);

    protected (IEnumerable<Point2D>, FiniteElement[]) BaseBuild(MeshParameters parameters)
    {
        var result = new
        {
            Points = new Point2D[(parameters.SplitsX.Sum() + 1) * (parameters.SplitsY.Sum() + 1)],
            Elements = new FiniteElement[parameters.SplitsX.Sum() * parameters.SplitsY.Sum()]
        };

        double[] pointsX = new double[parameters.SplitsX.Sum() + 1];
        double[] pointsY = new double[parameters.SplitsY.Sum() + 1];

        pointsX[0] = parameters.LinesX[0];
        pointsY[0] = parameters.LinesY[0];

        var idx = 1;

        for (int i = 0; i < parameters.LinesX.Length - 1; i++)
        {
            var hx = (parameters.LinesX[i + 1] - parameters.LinesX[i]) / parameters.SplitsX[i];

            for (int j = 0, k = i * parameters.SplitsX[i] + 1; j < parameters.SplitsX[i]; j++, k++)
            {
                pointsX[idx++] = pointsX[k - 1] + hx;
            }
        }

        idx = 1;

        for (int i = 0; i < parameters.LinesY.Length - 1; i++)
        {
            var hy = (parameters.LinesY[i + 1] - parameters.LinesY[i]) / parameters.SplitsY[i];

            for (int j = 0, k = i * parameters.SplitsX[i] + 1; j < parameters.SplitsY[i]; j++, k++)
            {
                pointsY[idx++] = pointsY[k - 1] + hy;
            }
        }

        idx = 0;

        for (int j = 0; j < pointsX.Length; j++)
        {
            for (int i = 0; i < pointsY.Length; i++)
            {
                result.Points[idx++] = new(pointsX[i], pointsY[j]);
            }
        }

        int nx = pointsX.Length;
        idx = 0;

        var nodes = new int[ElementSize];

        for (int j = 0; j < pointsY.Length - 1; j++)
        {
            for (int i = 0; i < pointsX.Length - 1; i++)
            {
                nodes[0] = i + j * nx;
                nodes[1] = i + 1 + j * nx;
                nodes[2] = i + (j + 1) * nx;
                nodes[3] = i + 1 + (j + 1) * nx;

                result.Elements[idx++] = new(nodes.ToArray(), FindNumberArea(result.Points, nodes, parameters));
            }
        }

        return (result.Points, result.Elements);
    }

    protected static int FindNumberArea(Point2D[] points, IEnumerable<int> nodes, MeshParameters parameters)
    {
        var localPoints = nodes.Select(node => points[node]).ToArray();

        foreach (var area in from area in parameters.Areas
                 let massCenter = new Point2D(localPoints.Sum(p => p.X) / 4.0, localPoints.Sum(p => p.Y) / 4.0)
                 where massCenter.X >= parameters.LinesX[area.X1] && massCenter.X <= parameters.LinesX[area.X2] &&
                       massCenter.Y >= parameters.LinesY[area.Y1] && massCenter.Y <= parameters.LinesY[area.Y2]
                 select area)
        {
            return area.Number;
        }

        throw new("Incorrect area parameters!");
    }
}

public class LinearMeshBuilder : MeshBuilder
{
    protected override int ElementSize => 4;

    public override (IEnumerable<Point2D>, FiniteElement[]) Build(MeshParameters parameters) => BaseBuild(parameters);
}

public class SuperMesh : IBaseMesh
{
    private readonly IEnumerable<Point2D> _points;
    private readonly FiniteElement[] _elements;

    public IReadOnlyList<Point2D> Points => _points.ToList().AsReadOnly();
    public IReadOnlyList<FiniteElement> Elements => _elements;
    public IReadOnlyList<Area> Areas { get; }

    public SuperMesh(MeshParameters parameters, MeshBuilder meshBuilder)
        => ((_points, _elements), Areas) = (meshBuilder.Build(parameters), parameters.Areas);
}

public readonly record struct Area(int Number, double Permeability, double Current, int X1, int X2, int Y1, int Y2)
{
    private const double VacuumPermeability = 4.0 * Math.PI * 10E-07;

    public static Area Parse(string line)
    {
        if (!TryParse(line, out var area))
        {
            throw new FormatException("Cant parse Area!");
        }

        return area;
    }

    public static bool TryParse(string line, out Area area)
    {
        var data = line.Split(new[] { ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);

        if (data.Length != 7 || !int.TryParse(data[0], out var number) || !double.TryParse(data[1], out var mu)
            || !double.TryParse(data[2], out var current) ||
            !int.TryParse(data[3], out var x1) || !int.TryParse(data[4], out var x2) ||
            !int.TryParse(data[5], out var y1) || !int.TryParse(data[6], out var y2))
        {
            area = default;
            return false;
        }

        area = new(number, VacuumPermeability * mu, current, x1, x2, y1, y2);
        return true;
    }
}

public class MeshParameters
{
    public ImmutableArray<double> LinesX { get; init; }
    public ImmutableArray<double> LinesY { get; init; }
    public ImmutableArray<int> SplitsX { get; init; }
    public ImmutableArray<int> SplitsY { get; init; }
    public ImmutableArray<int> Kx { get; init; }
    public ImmutableArray<int> Ky { get; init; }
    public (int, int) Nesting { get; init; }
    public ImmutableArray<Area> Areas { get; init; }

    public MeshParameters(IEnumerable<double> linesX, IEnumerable<double> linesY, IEnumerable<int> splitsX,
        IEnumerable<int> splitsY, IEnumerable<int> kx, IEnumerable<int> ky, (int, int) nesting,
        IEnumerable<Area> areas)
    {
        LinesX = linesX.ToImmutableArray();
        LinesY = linesY.ToImmutableArray();
        SplitsX = splitsX.ToImmutableArray();
        SplitsY = splitsY.ToImmutableArray();
        Kx = kx.ToImmutableArray();
        Ky = ky.ToImmutableArray();
        Nesting = nesting;
        Areas = areas.ToImmutableArray();
    }


    public MeshParameters(string path)
    {
        if (!File.Exists(path))
        {
            throw new("File does not exist");
        }

        using var sr = new StreamReader(path);
        LinesX = sr.ReadLine()!.Split().Select(double.Parse).ToImmutableArray();
        LinesY = sr.ReadLine()!.Split().Select(double.Parse).ToImmutableArray();
        SplitsX = sr.ReadLine()!.Split().Select(int.Parse).ToImmutableArray();
        SplitsY = sr.ReadLine()!.Split().Select(int.Parse).ToImmutableArray();
        Kx = sr.ReadLine()!.Split().Select(int.Parse).ToImmutableArray();
        Ky = sr.ReadLine()!.Split().Select(int.Parse).ToImmutableArray();
        var line = sr.ReadLine()!.Split();
        Nesting = (int.Parse(line[0]), int.Parse(line[1]));
        Areas = sr.ReadToEnd().Split("\n").Select(Area.Parse).ToImmutableArray();
    }
}